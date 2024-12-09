use gpf::{boolean3d::{boolean3d, BooleanType, SimpleBody}, geometry::{BBox, Cylinder, Plane, Sphere, Surf}};

#[test]
fn test1() {
	let body1 = {
		let left = Surf::Plane(Plane::new(-0.5, 0.0, 0.0, -1.0, 0.0, 0.0));
		let right = Surf::Cylinder(Cylinder::new(-0.5, 0.0, 0.0, 0.9, 0.3, (1.0f64 - 0.81 - 0.09).sqrt(), 1.1));
		let bottom = Surf::Plane(Plane::new(0.0, 0.0, -0.5, 0.0, 0.0, -1.0));
		let top = Surf::Cylinder(Cylinder::new(0.0, 0.0, -0.5, 0.0, 1.0, 0.0, 1.0));
		let front = Surf::Plane(Plane::new(0.0, -0.5, 0.0, 0.0, -1.0, 0.0));
		let back = Surf::Plane(Plane::new(0.0, 0.5, 0.0, 0.0, 1.0, 0.0));
		let bbox = BBox::new(-0.5, -0.5, -0.5, 1.0, 1.0, 1.0);
		SimpleBody::new(vec![left, right, bottom, top, front, back], bbox)
	};

	let body2 = {
		let left = Surf::Plane(Plane::new(0.0, 0.0, 0.0, -1.0, 0.0, 0.0));
		let right = Surf::Plane(Plane::new(1.0, 0.0, 0.0, 1.0, 0.0, 0.0));
		let bottom = Surf::Plane(Plane::new(0.0, 0.0, 0.0, 0.0, 0.0, -1.0));
		let top = Surf::Plane(Plane::new(0.0, 0.0, 1.0, 0.0, 0.0, 1.0));
		let front = Surf::Plane(Plane::new(0.0, 0.0, 0.0, 0.0, -1.0, 0.0));
		let back = Surf::Plane(Plane::new(0.0, 1.0, 0.0, 0.0, 1.0, 0.0));
		let bbox = BBox::new(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
		SimpleBody::new(vec![left, right, bottom, top, front, back], bbox)
	};
	boolean3d(&body1, &body2, BooleanType::Union, 0.001);
}

#[test]
fn test_component() {
    let bbox = BBox::new(-0.2, -0.2, -0.2, 1.0, 1.0, 1.0);
    let body1 = SimpleBody::new(vec![
        Surf::Sphere(Sphere::new(0.0, 0.0, 0.0, 0.2)),
        Surf::Sphere(Sphere::new(0.0, 0.0, 0.0, 0.1)),
        Surf::Sphere(Sphere::new(0.5, 0.0, 0.0, 0.2)),
        Surf::Sphere(Sphere::new(0.4, 0.0, 0.0, 0.2)),
    ], bbox);
    let body2 = SimpleBody::new(vec![], BBox::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

    boolean3d(&body1, &body2, BooleanType::Union, 0.0001);
}
