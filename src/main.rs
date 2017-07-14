extern crate argparse;
extern crate byteorder;

use byteorder::{LittleEndian, WriteBytesExt};
use field::{BoundaryCondition, Dimensions, Field, GridPoint};
use output::{FileWriter, OutputFormat};
use std::f64;
use std::fs::File;
use std::io::Write;
use std::ops::Mul;

mod output;
mod field;
mod parse;

pub type Vector = (f64, f64, f64);

fn gauss_3d(xx: f64, yy: f64, zz: f64, xm: Vector, dxm: Vector, f: f64) -> f64 {
    let r = (xx - xm.0, yy - xm.1, zz - xm.2);
    let exponent = -sqr(r.0 / dxm.0) - sqr(r.1 / dxm.1) - sqr(r.2 / dxm.2);
    exponent.exp() * (f * 2.0 * f64::consts::PI * r.0).sin()
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

const SIGMA_PML: Vector = (30.0, 30.0, 30.0);
const OUTPUT_EVERY: i32 = 10;
const GHOST: Dimensions = (1, 1, 1);
const LENGTH: Vector = (1.0, 1.0, 1.0);

fn main() {
    let parameters = parse::parse_parameters();

    if parameters.verbose {
        println!("{:#?}", parameters);
    }

    let grid = (parameters.grid, parameters.grid, parameters.grid);
    let sigma_j = (parameters.conductivity, parameters.conductivity, parameters.conductivity);

    let mut file = File::create(parameters.output_name).unwrap();

    let grid_0 = match parameters.output_format {
        OutputFormat::Bin2d |
        OutputFormat::Asc2d => 1,
        OutputFormat::Bin3d |
        OutputFormat::Asc3d => grid.0,
    };

    let to_write = [parameters.steps, grid_0, grid.1, grid.2, OUTPUT_EVERY];

    match parameters.output_format {
        OutputFormat::Asc3d |
        OutputFormat::Asc2d => {
            for element in &to_write {
                writeln!(file, "{}", element).unwrap();
            }
        }
        OutputFormat::Bin2d |
        OutputFormat::Bin3d => {
            for element in &to_write {
                file.write_i32::<LittleEndian>(*element as i32).unwrap();
            }
        }	
    };

    let dx = (LENGTH.0 / grid.0 as f64, LENGTH.1 / grid.1 as f64, LENGTH.2 / grid.2 as f64);
    let inv_sum = 1.0 / sqr(dx.0) + 1.0 / sqr(dx.1) + 1.0 / sqr(dx.2);
    let dt = 0.75 / inv_sum.sqrt();

    let x0 = (0.2, 0.8, 0.5);
    let sigma_gauss = (0.1, 0.1, 0.1);
    let freq = 10.0;

    let mut field = Field::create(grid, GHOST);
    field::foreach_3d((-GHOST.0, -GHOST.1, -GHOST.2), (grid.0 + GHOST.0, grid.1 + GHOST.1, grid.2 + GHOST.2), |ix, iy, iz| {
        field[(ix, iy, iz)].e.2 = gauss_3d(ix as f64 * dx.0, iy as f64 * dx.1, (iz as f64 + 0.5) * dx.2, x0, sigma_gauss, freq);
        field[(ix, iy, iz)].h.1 = -gauss_3d((ix as f64 + 0.5) * dx.0, iy as f64 * dx.1, (iz as f64 + 0.5) * dx.2, x0, sigma_gauss, freq);
    });

    let cn = (dt / dx.0, dt / dx.1, dt / dx.2);

    /*
	init_profs(fdtd_integration_loop);
	init_profs(e_field_update);
	init_profs(h_field_update);
	init_profs(output);

	profs_start(fdtd_integration_loop);
	*/

    let mut file_writer = FileWriter::new(file, parameters.output_format);

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
                SIGMA_PML
            } else {
                sigma_j
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
                SIGMA_PML
            } else {
                sigma_j
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
            file_writer.write_field(&mut field);
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
