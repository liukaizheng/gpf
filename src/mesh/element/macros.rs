#[macro_export]
macro_rules! element_id {
    (struct $name: ident) => {
        impl From<usize> for $name {
            #[inline(always)]
            fn from(id: usize) -> Self {
                Self(id)
            }
        }

        impl Default for $name {
            #[inline(always)]
            fn default() -> Self {
                Self(INVALID_IND)
            }
        }

        impl<T> Index<$name> for Vec<T> {
            type Output = T;
            #[inline(always)]
            fn index(&self, index: $name) -> &Self::Output {
                &self[index.0]
            }
        }

        impl<T> IndexMut<$name> for Vec<T> {
            #[inline(always)]
            fn index_mut(&mut self, index: $name) -> &mut Self::Output {
                &mut self[index.0]
            }
        }

        impl<'b, T> Index<$name> for bumpalo::collections::Vec<'b, T> {
            type Output = T;
            #[inline(always)]
            fn index(&self, index: $name) -> &Self::Output {
                &self[index.0]
            }
        }

        impl<'b, T> IndexMut<$name> for bumpalo::collections::Vec<'b, T> {
            #[inline(always)]
            fn index_mut(&mut self, index: $name) -> &mut Self::Output {
                &mut self[index.0]
            }
        }
        impl<T> Index<$name> for [T] {
            type Output = T;
            #[inline(always)]
            fn index(&self, index: $name) -> &Self::Output {
                &self[index.0]
            }
        }

        impl<T> IndexMut<$name> for [T] {
            #[inline(always)]
            fn index_mut(&mut self, index: $name) -> &mut Self::Output {
                &mut self[index.0]
            }
        }

        impl ElementIndex for $name {
            #[inline(always)]
            fn index(&self) -> usize {
                self.0
            }
        }
        impl ElementId for $name {}

        impl Deref for $name {
            type Target = usize;
            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }

        impl DerefMut for $name {
            fn deref_mut(&mut self) -> &mut Self::Target {
                &mut self.0
            }
        }
    };
}
