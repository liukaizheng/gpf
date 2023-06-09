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
            #[inline(always)]
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

            #[inline(always)]
            fn next(&mut self) -> Option<Self::Item> {
                iter_next(self)
            }
        }

    }
}

#[macro_export]
macro_rules! halfedges_iterator {
    (struct $name:ident) => {
        element_iterator! {
            struct $name -> HalfedgeId, {
                fn id(&self) -> Self::Id {
                    self.curr_he
                }
            }, {
                fn valid(&self) -> bool {
                    true
                }
            }, {
                fn next(&mut self) {
                    self.just_start = false;
                    self.next();
                }
            }, {
                fn is_end(&self) -> bool {
                    !self.just_start && self.curr_he == self.first_he
                }
            }

        }
    };
}
