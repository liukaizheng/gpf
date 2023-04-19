#[macro_export]
macro_rules! element_iterator {
    (struct $name: ident -> $item: ty, {$($id: tt)*}, {$($valid: tt)*}, {$($next: tt)*}) => {
        impl<'a, M: Mesh> Element for $name<'a, M> {
            type Item = $item;
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