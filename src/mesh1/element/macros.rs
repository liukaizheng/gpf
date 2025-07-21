#[macro_export]
macro_rules! element_iter_struct {
    (struct $name:ident -> $mesh_: tt, $elem_id:ty, $elem:tt, $from_ref:ident, $into_ref:ident, $elem_data: ident, {$( $mut_:tt )?}) => {
        pub struct $name<'m, M: $mesh_> {
            pub id: $elem_id,
            pub data: &'m $($mut_)? M::$elem,
            pub mesh: NonNull<M>,
            _marker: PhantomData<&'m $($mut_)? M>,
        }

        impl<'m, M: $mesh_> $name<'m, M> {
            pub fn new(id: $elem_id, mesh: &'m $($mut_)? M) -> Self {
                let $($mut_)? mesh = NonNull::$from_ref(mesh);
                unsafe {
                    let data = mesh.$into_ref().$elem_data(id);
                    Self {
                        id,
                        data,
                        mesh,
                        _marker: PhantomData,
                    }
                }
            }

            pub fn new_with_data(id: $elem_id, data: &'m $($mut_)? M::$elem, mesh: &'m $($mut_)? M) -> Self {
                let $($mut_)? mesh = NonNull::$from_ref(mesh);
                unsafe {
                    Self {
                        id,
                        data,
                        mesh,
                        _marker: PhantomData,
                    }
                }
            }
        }
    };

}
