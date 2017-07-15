use argparse::{ArgumentParser, Store, StoreConst, StoreTrue};
use field::BoundaryCondition;
use output::OutputFormat;
use std::cmp;

#[derive(Debug)]
pub struct Parameters {
    pub grid: i32,
    pub steps: i32,
    pub conductivity: f64,
    pub boundary_condition: BoundaryCondition,
    pub output_format: OutputFormat,
    pub output_name: String,
    pub verbose: bool,
}

pub fn parse_parameters() -> Parameters {
    let mut parameters = Parameters {
        grid: 50,
        steps: 100,
        conductivity: 0.0,
        boundary_condition: BoundaryCondition::Reflecting,
        output_format: OutputFormat::Bin3d,
        output_name: String::from("out"),
        verbose: false,
    };

    let mut output_format = String::new();

    {
        let mut parser = ArgumentParser::new();
        parser.refer(&mut parameters.grid).add_option(
            &["-g", "--grid"],
            Store,
            "set grid cell number",
        );
        parser.refer(&mut parameters.steps).add_option(
            &["-t", "--steps"],
            Store,
            "set number of timesteps",
        );
        parser.refer(&mut parameters.conductivity).add_option(
            &[
                "-s",
                "--conductivity",
            ],
            Store,
            "set conductivity in box",
        );
        parser
            .refer(&mut parameters.boundary_condition)
            .add_option(&["-p", "--periodic"], StoreConst(BoundaryCondition::Periodic), "set periodic boundaries")
            .add_option(&["-r", "--reflecting"], StoreConst(BoundaryCondition::Reflecting), "set reflecting boundaries")
            .add_option(&["-a", "--absorbing"], StoreConst(BoundaryCondition::Absorbing), "set absorbing boundaries (pmls)");
        parser.refer(&mut output_format).add_option(
            &["-o", "--output"],
            Store,
            "set type of output",
        );
        parser.refer(&mut parameters.output_name).add_option(
            &["-n", "--name"],
            Store,
            "set name for outputfile",
        );
        parser.refer(&mut parameters.verbose).add_option(
            &["-v", "--verbose"],
            StoreTrue,
            "set to verbose (profiling output)",
        );
        parser.parse_args_or_exit();
    }

    parameters.grid = cmp::max(parameters.grid, 1);
    parameters.steps = cmp::max(parameters.steps, 1);
    parameters.conductivity = f64::max(parameters.conductivity, 0.0);

    match OutputFormat::parse(&output_format) {
        Ok(output_format) => parameters.output_format = output_format,
        Err(message) => println!("{}\nFalling back to {:?}", message, parameters.output_format),
    }

    parameters
}