use crate::field::{BoundaryCondition, Dimensions, Field, GridPoint};
use crate::output::{OutputFormat, OutputWriter};
use std::f64;
use std::fs::File;
use std::io::BufWriter;
use std::ops::Mul;

mod output;
mod field;
mod parse;

pub type Vector = (f64, f64, f64);

fn gauss_3d(x: Vector, xm: Vector, dxm: Vector, frequency_in_x_direction: f64) -> f64 {
    let r = (x.0 - xm.0, x.1 - xm.1, x.2 - xm.2);
    let exponent = -sqr(r.0 / dxm.0) - sqr(r.1 / dxm.1) - sqr(r.2 / dxm.2);
    exponent.exp() * (2.0 * f64::consts::PI * frequency_in_x_direction * r.0).sin()
}

fn sqr<T: Copy + Mul<T>>(value: T) -> T::Output {
    value * value
}

fn exchange_hfields(field: &mut Field, grid: Dimensions) {
    for z in 0..grid.2 {
        for y in 0..grid.1 {
            field[(-1, y, z)].h = field[(grid.0 - 1, y, z)].h;
            field[(grid.0, y, z)].h = field[(0, y, z)].h;
        }
    }

    for z in 0..grid.2 {
        for x in 0..grid.0 {
            field[(x, -1, z)].h = field[(x, grid.1 - 1, z)].h;
            field[(x, grid.1, z)].h = field[(x, 0, z)].h;
        }
    }

    for y in 0..grid.1 {
        for x in 0..grid.0 {
            field[(x, y, -1)].h = field[(x, y, grid.2 - 1)].h;
            field[(x, y, grid.2)].h = field[(x, y, 0)].h;
        }
    }
}

fn exchange_efields(field: &mut Field, grid: Dimensions) {
    for z in 0..grid.2 {
        for y in 0..grid.1 {
            field[(-1, y, z)].e = field[(grid.0 - 1, y, z)].e;
            field[(grid.0, y, z)].e = field[(0, y, z)].e;
        }
    }

    for z in 0..grid.2 {
        for x in 0..grid.0 {
            field[(x, -1, z)].e = field[(x, grid.1 - 1, z)].e;
            field[(x, grid.1, z)].e = field[(x, 0, z)].e;
        }
    }

    for y in 0..grid.1 {
        for x in 0..grid.0 {
            field[(x, y, -1)].e = field[(x, y, grid.2 - 1)].e;
            field[(x, y, grid.2)].e = field[(x, y, 0)].e;
        }
    }
}

const PML_CONDUCTIVITY: Vector = (30.0, 30.0, 30.0);
const OUTPUT_EVERY: i32 = 10;
const GHOST: Dimensions = (1, 1, 1);
const LENGTH: Vector = (1.0, 1.0, 1.0);

fn main() {
    let parameters = parse::parse_parameters();

    if parameters.verbose {
        println!("{:#?}", parameters);
    }

    let grid = (parameters.grid, parameters.grid, parameters.grid);
    let conductivity = (parameters.conductivity, parameters.conductivity, parameters.conductivity);

    let writer = match parameters.output_format {
        OutputFormat::Asc3d |
        OutputFormat::Asc2d => OutputWriter::Ascii,
        OutputFormat::Bin2d |
        OutputFormat::Bin3d => OutputWriter::Binary,
    };

    let grid_0 = match parameters.output_format {
        OutputFormat::Bin2d |
        OutputFormat::Asc2d => 1,
        OutputFormat::Bin3d |
        OutputFormat::Asc3d => grid.0,
    };

    let file = File::create(parameters.output_name).unwrap();
    let mut target = BufWriter::new(file);
    writer.write_metadata(&mut target, &[parameters.steps, grid_0, grid.1, grid.2, OUTPUT_EVERY]);

    let dx = (LENGTH.0 / f64::from(grid.0), LENGTH.1 / f64::from(grid.1), LENGTH.2 / f64::from(grid.2));
    let inv_sum = 1.0 / sqr(dx.0) + 1.0 / sqr(dx.1) + 1.0 / sqr(dx.2);
    let dt = 0.75 / inv_sum.sqrt();

    let x0 = (0.2, 0.8, 0.5);
    let sigma_gauss = (0.1, 0.1, 0.1);
    let freq = 10.0;

    let mut field = Field::create(grid, GHOST);
    field::foreach_3d((-GHOST.0, -GHOST.1, -GHOST.2), (grid.0 + GHOST.0, grid.1 + GHOST.1, grid.2 + GHOST.2), |ix, iy, iz| {
        field[(ix, iy, iz)].e.2 = gauss_3d((f64::from(ix) * dx.0, f64::from(iy) * dx.1, (f64::from(iz) + 0.5) * dx.2), x0, sigma_gauss, freq);
        field[(ix, iy, iz)].h.1 = -gauss_3d(((f64::from(ix) + 0.5) * dx.0, f64::from(iy) * dx.1, (f64::from(iz) + 0.5) * dx.2), x0, sigma_gauss, freq);
    });

    let cn = (dt / dx.0, dt / dx.1, dt / dx.2);

    /*
	init_profs(fdtd_integration_loop);
	init_profs(e_field_update);
	init_profs(h_field_update);
	init_profs(output);

	profs_start(fdtd_integration_loop);
	*/

    let boundary_condition = parameters.boundary_condition;
    for t in 0..parameters.steps {
        if boundary_condition == BoundaryCondition::Reflecting {
            exchange_hfields(&mut field, grid);
        }

        /*
		if ((verb) && ((t%output_every) == 0)) profs_start(e_field_update);
		*/

        field::foreach_3d((0, 0, 0), (grid.0 + 1, grid.1 + 1, grid.2 + 1), |ix, iy, iz| {
            let sigma = if boundary_condition == BoundaryCondition::Absorbing && is_inside_pml((ix, iy, iz), grid) {
                PML_CONDUCTIVITY
            } else {
                conductivity
            };
            field[(ix, iy, iz)].e = (
                (cn.1 * (field[(ix, iy, iz)].h.2 - field[(ix, iy - 1, iz)].h.2) - cn.2 * (field[(ix, iy, iz)].h.1 - field[(ix, iy, iz - 1)].h.1) +
                     field[(ix, iy, iz)].e.0 * (1. - sigma.0 / 2.0 * dt)) / (1. + sigma.0 / 2.0 * dt),
                (cn.2 * (field[(ix, iy, iz)].h.0 - field[(ix, iy, iz - 1)].h.0) - cn.0 * (field[(ix, iy, iz)].h.2 - field[(ix - 1, iy, iz)].h.2) +
                     field[(ix, iy, iz)].e.1 * (1. - sigma.1 / 2.0 * dt)) / (1. + sigma.1 / 2.0 * dt),
                (cn.0 * (field[(ix, iy, iz)].h.1 - field[(ix - 1, iy, iz)].h.1) - cn.1 * (field[(ix, iy, iz)].h.0 - field[(ix, iy - 1, iz)].h.0) +
                     field[(ix, iy, iz)].e.2 * (1. - sigma.2 / 2.0 * dt)) / (1. + sigma.2 / 2.0 * dt),
            );
        });

        /*
		if ((verb) && ((t%output_every) == 0)) { profs_end(e_field_update);}
		*/

        if boundary_condition == BoundaryCondition::Reflecting {
            exchange_efields(&mut field, grid);
        }

        /*
		if ((verb) && ((t%output_every) == 0)) profs_start(h_field_update);
		*/

        field::foreach_3d((-1, -1, -1), grid, |ix, iy, iz| {
            let sigma = if boundary_condition == BoundaryCondition::Absorbing && is_inside_pml((ix, iy, iz), grid) {
                PML_CONDUCTIVITY
            } else {
                conductivity
            };
            field[(ix, iy, iz)].h = (
                (cn.1 * (field[(ix, iy, iz)].e.2 - field[(ix, iy + 1, iz)].e.2) - cn.2 * (field[(ix, iy, iz)].e.1 - field[(ix, iy, iz + 1)].e.1) +
                     field[(ix, iy, iz)].h.0 * (1. - sigma.0 / 2.0 * dt)) / (1. + sigma.0 / 2.0 * dt),
                (cn.2 * (field[(ix, iy, iz)].e.0 - field[(ix, iy, iz + 1)].e.0) - cn.0 * (field[(ix, iy, iz)].e.2 - field[(ix + 1, iy, iz)].e.2) +
                     field[(ix, iy, iz)].h.1 * (1. - sigma.1 / 2.0 * dt)) / (1. + sigma.1 / 2.0 * dt),
                (cn.0 * (field[(ix, iy, iz)].e.1 - field[(ix + 1, iy, iz)].e.1) - cn.1 * (field[(ix, iy, iz)].e.0 - field[(ix, iy + 1, iz)].e.0) +
                     field[(ix, iy, iz)].h.2 * (1. - sigma.2 / 2.0 * dt)) / (1. + sigma.2 / 2.0 * dt),
            )
        });

        /*
		if ((verb) && ((t%output_every) == 0)) { profs_end(h_field_update);}

		if ((verb) && ((t%output_every) == 0)) {profs_start(output);}
		*/

        if t % OUTPUT_EVERY == 0 {
            writer.write_field(&mut target, &field);
        }

        /*
		if ((verb) && ((t%output_every) == 0)) { profs_end(output); }
		*/
    }

    /*
	profs_end(fdtd_integration_loop);
	*/
}

const PML: Dimensions = (7, 7, 7);

fn is_inside_pml(point: GridPoint, grid: Dimensions) -> bool {
    point.0 < PML.0 || point.0 >= grid.0 - PML.0 || point.1 < PML.1 || point.1 >= grid.1 - PML.1 || point.2 < PML.2 || point.2 >= grid.2 - PML.2
}
