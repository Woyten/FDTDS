use byteorder::{LittleEndian, WriteBytesExt};
use field::{self, Field};
use std::fs::File;
use std::io::BufWriter;

#[derive(Debug, Eq, PartialEq)]
pub enum OutputFormat {
    Asc2d,
    Asc3d,
    Bin2d,
    Bin3d,
}

impl OutputFormat {
    pub fn parse(key: &str) -> Result<OutputFormat, &'static str> {
        match key {
            "BIN3D" | "bin3d" | "binary3d" => Ok(OutputFormat::Bin3d),
            "BIN2D" | "bin2d" | "binary2d" => Ok(OutputFormat::Bin2d),
            "ASC3D" | "asc3d" | "ascii3d" => Ok(OutputFormat::Asc3d),
            "ASC2D" | "asc2d" | "ascii2d" => Ok(OutputFormat::Asc2d),
            _ => Err("no valid output type specified: options are ASC/asc/ascii for textual output and BIN/bin/binary for binary"),
        }
    }
}

/*
void write_text_3d(struct field* fp,int* im, FILE* file, int output_every, int t, int e){
                if ((t%output_every) == 0) {
                        foreach_3d(ix, iy, iz, 0, 0) {
                                fprintf(file,"%.9g\n",F3(fp, e, ix, iy, iz));
                        } foreach_3d_end;
                }
}

void write_binary_2d(struct field* fp,int* im, FILE* file, int output_every, int t, int e){
                if ((t%output_every) == 0) {
                        foreach_2d(ix, iy, 0, 0){
                                fwrite(&F3(fp, e, ix, iy, im[2]/2), sizeof(double), 1, file);
                        } foreach_2d_end;
                }
}

void write_text_2d(struct field* fp,int* im, FILE* file, int output_every, int t, int e){
                if ((t%output_every) == 0) {
                        foreach_2d(ix, iy, 0, 0){
                                fprintf(file,"%.9g\n",F3(fp, e, ix, iy, im[2]/2));
                        } foreach_2d_end;
                }
}
 */

pub struct FileWriter {
    file: BufWriter<File>,
    output_format: OutputFormat,
}

impl FileWriter {
    pub fn new(file: File, output_format: OutputFormat) -> Self {
        FileWriter {
            file: BufWriter::new(file),
            output_format,
        }
    }

    pub fn write_field(&mut self, field: &mut Field) {
        let ghost = field.ghost;
        let total_grid = field.total_grid;
        field::foreach_3d((0, 0, 0), (total_grid.0 - 2 * ghost.0, total_grid.1 - 2 * ghost.1, total_grid.2 - 2 * ghost.2), |ix, iy, iz| {
            self.file
                .write_f64::<LittleEndian>(field.get(ix, iy, iz).e.2)
                .unwrap();
        });
    }
}