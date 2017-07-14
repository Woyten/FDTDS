use Vector;
use std::ops::{Index, IndexMut};

pub type Dimensions = (i32, i32, i32);
pub type GridPoint = (i32, i32, i32);

#[derive(Clone, Default)]
pub struct FieldStrength {
    pub e: Vector,
    pub h: Vector,
}

pub struct Field {
    data: Vec<FieldStrength>,
    pub total_grid: Dimensions,
    pub ghost: Dimensions,
}

impl Field {
    pub fn create(grid: Dimensions, ghost: Dimensions) -> Self {
        let total_grid = (grid.0 + 2 * ghost.0, grid.1 + 2 * ghost.1, grid.2 + 2 * ghost.2);

        Field {
            data: vec![FieldStrength::default(); total_grid.0 as usize * total_grid.1 as usize * total_grid.2 as usize],
            total_grid,
            ghost,
        }
    }

    fn get_index_for_point(&self, grid_point: GridPoint) -> usize {
        let z = grid_point.2 + self.ghost.2;
        let yz = grid_point.1 + self.ghost.1 + z * self.total_grid.1;
        let xyz = grid_point.0 + self.ghost.0 + yz * self.total_grid.0;
        xyz as usize
    }
}

impl Index<GridPoint> for Field {
    type Output = FieldStrength;

    fn index(&self, grid_point: GridPoint) -> &Self::Output {
        let index = self.get_index_for_point(grid_point);
        unsafe { self.data.get_unchecked(index) }
    }
}

impl IndexMut<GridPoint> for Field {
    fn index_mut(&mut self, grid_point: GridPoint) -> &mut Self::Output {
        let index = self.get_index_for_point(grid_point);
        unsafe { self.data.get_unchecked_mut(index) }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
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