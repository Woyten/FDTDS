use byteorder::{LittleEndian, WriteBytesExt};
use field::{self, Field};
use std::io::Write;

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

pub enum OutputWriter {
    Binary,
    Ascii,
}

impl OutputWriter {
    pub fn write_metadata<W: Write>(&self, target: &mut W, values: &[i32]) {
        match *self {
            OutputWriter::Binary => {
                for value in values {
                    target.write_i32::<LittleEndian>(*value).unwrap();
                }
            }
            OutputWriter::Ascii => {
                for value in values {
                    writeln!(target, "{}", value).unwrap();
                }
            }
        }
    }

    pub fn write_field<W: Write>(&self, target: &mut W, field: &Field) {
        let ghost = field.ghost;
        let total_grid = field.total_grid;
        match *self {
            OutputWriter::Binary => {
                field::foreach_3d((0, 0, 0), (total_grid.0 - 2 * ghost.0, total_grid.1 - 2 * ghost.1, total_grid.2 - 2 * ghost.2), |ix, iy, iz| {
                    target
                        .write_f64::<LittleEndian>(field[(ix, iy, iz)].e.2)
                        .unwrap();
                });
            }
            OutputWriter::Ascii => {
                field::foreach_3d((0, 0, 0), (total_grid.0 - 2 * ghost.0, total_grid.1 - 2 * ghost.1, total_grid.2 - 2 * ghost.2), |ix, iy, iz| {
                    writeln!(target, "{:.8e}", field[(ix, iy, iz)].e.2).unwrap();
                });
            }
        }
    }
}