#[macro_export]
macro_rules! element_id {
    (struct $name: ident) => {
        impl From<usize> for $name {
            #[inline(always)]
            fn from(id: usize) -> Self {
                Self(id)
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

        impl ElementId for $name {
            #[inline(always)]
            fn new() -> Self {
                Self(INVALID_IND)
            }
            #[inline(always)]
            fn valid(&self) -> bool {
                self.0 != INVALID_IND
            }
        }

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
#[macro_export]
macro_rules! element_iterator {
    (struct $name: ident -> $item: ty, {$($id: tt)*}, {$($valid: tt)*}, {$($next: tt)*}, {$($is_end: tt)*}) => {
        impl<'a, M: Mesh> Element for $name<'a, M> {
            type Id = $item;
            type M = M;
            #[inline]
            $($id)*

            #[inline(always)]
            fn mesh(&self) -> &M {
                self.mesh
            }

            #[inline(always)]
            $($valid)*

            #[inline]
            $($next)*

            #[inline(always)]
            $($is_end)*
        }
        impl<'a, M: Mesh> Iterator for $name<'a, M> {
            type Item = $item;

            #[inline]
            fn next(&mut self) -> Option<Self::Item> {
                iter_next(self)
            }
        }

    }
}
