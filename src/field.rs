pub type Dimensions = (i32, i32, i32);
pub type GridPoint = (i32, i32, i32);

pub struct Field {
    data: Vec<f64>,
    pub total_grid: Dimensions,
    pub ghost: Dimensions,
}

impl Field {
    pub fn create(grid: Dimensions, ghost: Dimensions) -> Self {
        let nr_components = 6;
        let total_grid = (grid.0 + 2 * ghost.0, grid.1 + 2 * ghost.1, grid.2 + 2 * ghost.2);

        Field {
            data: vec![0.0; nr_components as usize * total_grid.0 as usize * total_grid.1 as usize * total_grid.2 as usize],
            total_grid,
            ghost,
        }
    }

    pub fn e_x(&mut self, jx: i32, jy: i32, jz: i32) -> &mut f64 {
        let index = self.f3_off(0, jx, jy, jz);
        &mut self.data[index]
    }

    pub fn e_y(&mut self, jx: i32, jy: i32, jz: i32) -> &mut f64 {
        let index = self.f3_off(1, jx, jy, jz);
        &mut self.data[index]
    }

    pub fn e_z(&mut self, jx: i32, jy: i32, jz: i32) -> &mut f64 {
        let index = self.f3_off(2, jx, jy, jz);
        &mut self.data[index]
    }

    pub fn h_x(&mut self, jx: i32, jy: i32, jz: i32) -> &mut f64 {
        let index = self.f3_off(3, jx, jy, jz);
        &mut self.data[index]
    }

    pub fn h_y(&mut self, jx: i32, jy: i32, jz: i32) -> &mut f64 {
        let index = self.f3_off(4, jx, jy, jz);
        &mut self.data[index]
    }

    pub fn h_z(&mut self, jx: i32, jy: i32, jz: i32) -> &mut f64 {
        let index = self.f3_off(5, jx, jy, jz);
        &mut self.data[index]
    }

    fn f3_off(&self, fldnr: i32, jx: i32, jy: i32, jz: i32) -> usize {
        let index = jx + self.ghost.0 + self.total_grid.0 * (jy + self.ghost.1 + self.total_grid.1 * (jz + self.ghost.2 + self.total_grid.2 * fldnr));
        index as usize
    }
}

#[derive(Clone, Copy, Eq, PartialEq)]
pub enum BoundaryCondition {
    Reflecting,
    Periodic,
    Absorbing,
}

pub fn foreach_3d<A>(ilo: Dimensions, ihi: Dimensions, mut action: A)
where
    A: FnMut(i32, i32, i32),
{
    for iz in ilo.2..ihi.2 {
        for iy in ilo.1..ihi.1 {
            for ix in ilo.0..ihi.0 {
                action(ix, iy, iz);
            }
        }
    }
}