

			//******** BuffA Code Begins ********

vec2 R;float N;// R is resolution, N is number of samples
vec4 D ( vec2 U ) {return texture(BuffD,U/R);}//Samples previous field state
float Xa (vec2 U0, vec2 U, vec2 U1, inout vec4 Q, in vec2 r) {
    vec2 V = U + r, u = D(V).xy,
         V0 = V - u,
         V1 = V + u;
    float P = D (V0).z, rr = length(r);
    Q.xy -= r*(P-Q.z)/rr/N;
    return (0.5*(length(V0-U0)-length(V1-U1))+P)/N;
    /* This is the interaction function
		U1 is my current position
		U is my previous position
		Q is my previous state
		r is the relativity between me and my neighbor
    	V1 is my neighbor's current position
		V is my neighbor's previous position
		P is my neighbors pressure
		Q.z is my previous pressure
		rr is the distance between me and my neighbor
		r*(P-Q.z)/rr/N is this interaction's contribution
			to the my velocity. It is the gradient of
			the pressure between me and my neighbor
		length(V0-U0)-length(V1-U1) is the space contraction between
			the next and the previous state. Because pressure
			is energy per volume. If the volume contracts
			there must be a higher energy density. If
			the space expands, the energy is disipated
			over a larger area.
		I add P because I completely trade pressures 
			with my neighbors. This is why the best
			sampling patterns are like checkerboards.
			every other point in space should interact
			so that there is a feedback system.
			The physical explanation of this is that
			the fluid is kind of like many billiards
			balls. When two billiards balls collide,
			they completely swap energies.
		I divide the outputs by N to average each interaction.
    */
}

vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
   R = RENDERSIZE.xy;
 	vec2 U0 = U - D(U).xy, // get previous state
         U1 = U + D(U).xy; // get next state
 	float P = 0.; Q = D(U0);
 	N = 4.;// checkerboard sampling pattern
    P += Xa (U0,U,U1,Q, vec2( 1, 0) );
 	P += Xa (U0,U,U1,Q, vec2( 0,-1) );
 	P += Xa (U0,U,U1,Q, vec2(-1, 0) );
 	P += Xa (U0,U,U1,Q, vec2( 0, 1) );
 	Q.z = P;
 	
 	// Init and Walls
 	if (int(FRAMECOUNT) <= 1 || Reset != 0.) Q = vec4(0);
    if (U.x < 1.||U.y < 1.||R.x-U.x < 1.||R.y-U.y < 1.) Q.xy *= 0.;
    //This Jet setup will cause pressure to always increase in the system
 	if (length(U-vec2(0,0.5*R.y)) < 4.) {Q.xy= Q.xy*.9+.1*vec2(.5,0);; Q.w = 1.;}
    if (length(U-vec2(R.x,0.5*R.y)) < 4.) {Q.xy= Q.xy*.9+.1*vec2(-.5,0);; Q.w = 1.;}
	return Q; 
 } 


			//******** BuffB Code Begins ********

//vec2 R;float N;// R is resolution, N is number of samples
vec4 A ( vec2 U ) {return texture(BuffA,U/R);}//Samples previous field state
float Xb (vec2 U0, vec2 U, vec2 U1, inout vec4 Q, in vec2 r) {
    vec2 V = U + r, u = A(V).xy,
         V0 = V - u,
         V1 = V + u;
    float P = A (V0).z, rr = length(r);
    Q.xy -= r*(P-Q.z)/rr/N;
    return (0.5*(length(V0-U0)-length(V1-U1))+P)/N;
    /* This is the interaction function
		U1 is my current position
		U is my previous position
		Q is my previous state
		r is the relativity between me and my neighbor
    	V1 is my neighbor's current position
		V is my neighbor's previous position
		P is my neighbors pressure
		Q.z is my previous pressure
		rr is the distance between me and my neighbor
		r*(P-Q.z)/rr/N is this interaction's contribution
			to the my velocity. It is the gradient of
			the pressure between me and my neighbor
		length(V0-U0)-length(V1-U1) is the space contraction between
			the next and the previous state. Because pressure
			is energy per volume. If the volume contracts
			there must be a higher energy density. If
			the space expands, the energy is disipated
			over a larger area.
		I add P because I completely trade pressures 
			with my neighbors. This is why the best
			sampling patterns are like checkerboards.
			every other point in space should interact
			so that there is a feedback system.
			The physical explanation of this is that
			the fluid is kind of like many billiards
			balls. When two billiards balls collide,
			they completely swap energies.
		I divide the outputs by N to average each interaction.
    */
}

vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
   R = RENDERSIZE.xy;
 	vec2 U0 = U - A(U).xy, // get previous state
         U1 = U + A(U).xy; // get next state
 	float P = 0.; Q = A(U0);
 	N = 4.;// checkerboard sampling pattern
    P += Xb (U0,U,U1,Q, vec2( 1, 0) );
 	P += Xb (U0,U,U1,Q, vec2( 0,-1) );
 	P += Xb (U0,U,U1,Q, vec2(-1, 0) );
 	P += Xb (U0,U,U1,Q, vec2( 0, 1) );
 	Q.z = P;
 	
 	// Init and Walls
 	if (int(FRAMECOUNT) <= 1 || Reset != 0.) Q = vec4(0);
    if (U.x < 1.||U.y < 1.||R.x-U.x < 1.||R.y-U.y < 1.) Q.xy *= 0.;
    //This Jet setup will cause pressure to always increase in the system
 	if (length(U-vec2(0,0.5*R.y)) < 4.) {Q.xy= Q.xy*.9+.1*vec2(.5,0);; Q.w = 1.;}
    if (length(U-vec2(R.x,0.5*R.y)) < 4.) {Q.xy= Q.xy*.9+.1*vec2(-.5,0);; Q.w = 1.;}
	return Q; 
 } 


			//******** BuffC Code Begins ********

//vec2 R;float N;// R is resolution, N is number of samples
vec4 B ( vec2 U ) {return texture(BuffB,U/R);}//Samples previous field state
float Xc (vec2 U0, vec2 U, vec2 U1, inout vec4 Q, in vec2 r) {
    vec2 V = U + r, u = B(V).xy,
         V0 = V - u,
         V1 = V + u;
    float P = B (V0).z, rr = length(r);
    Q.xy -= r*(P-Q.z)/rr/N;
    return (0.5*(length(V0-U0)-length(V1-U1))+P)/N;
    /* This is the interaction function
		U1 is my current position
		U is my previous position
		Q is my previous state
		r is the relativity between me and my neighbor
    	V1 is my neighbor's current position
		V is my neighbor's previous position
		P is my neighbors pressure
		Q.z is my previous pressure
		rr is the distance between me and my neighbor
		r*(P-Q.z)/rr/N is this interaction's contribution
			to the my velocity. It is the gradient of
			the pressure between me and my neighbor
		length(V0-U0)-length(V1-U1) is the space contraction between
			the next and the previous state. Because pressure
			is energy per volume. If the volume contracts
			there must be a higher energy density. If
			the space expands, the energy is disipated
			over a larger area.
		I add P because I completely trade pressures 
			with my neighbors. This is why the best
			sampling patterns are like checkerboards.
			every other point in space should interact
			so that there is a feedback system.
			The physical explanation of this is that
			the fluid is kind of like many billiards
			balls. When two billiards balls collide,
			they completely swap energies.
		I divide the outputs by N to average each interaction.
    */
}

vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
   R = RENDERSIZE.xy;
 	vec2 U0 = U - B(U).xy, // get previous state
         U1 = U + B(U).xy; // get next state
 	float P = 0.; Q = B(U0);
 	N = 4.;// checkerboard sampling pattern
    P += Xc (U0,U,U1,Q, vec2( 1, 0) );
 	P += Xc (U0,U,U1,Q, vec2( 0,-1) );
 	P += Xc (U0,U,U1,Q, vec2(-1, 0) );
 	P += Xc (U0,U,U1,Q, vec2( 0, 1) );
 	Q.z = P;
 	
 	// Init and Walls
 	if (int(FRAMECOUNT) <= 1 || Reset != 0.) Q = vec4(0);
    if (U.x < 1.||U.y < 1.||R.x-U.x < 1.||R.y-U.y < 1.) Q.xy *= 0.;
    //This Jet setup will cause pressure to always increase in the system
 	if (length(U-vec2(0,0.5*R.y)) < 4.) {Q.xy= Q.xy*.9+.1*vec2(.5,0);; Q.w = 2.;}
    if (length(U-vec2(R.x,0.5*R.y)) < 4.) {Q.xy= Q.xy*.9+.1*vec2(-.5,0);; Q.w = 2.;}
	return Q; 
 } 


			//******** BuffD Code Begins ********

//vec2 R;float N;// R is resolution, N is number of samples
vec4 C ( vec2 U ) {return texture(BuffC,U/R);}//Samples previous field state
float Xd (vec2 U0, vec2 U, vec2 U1, inout vec4 Q, in vec2 r) {
    vec2 V = U + r, u = C(V).xy,
         V0 = V - u,
         V1 = V + u;
    float P = C (V0).z, rr = length(r);
    Q.xy -= r*(P-Q.z)/rr/N;
    return (0.5*(length(V0-U0)-length(V1-U1))+P)/N;
    /* This is the interaction function
		U1 is my current position
		U is my previous position
		Q is my previous state
		r is the relativity between me and my neighbor
    	V1 is my neighbor's current position
		V is my neighbor's previous position
		P is my neighbors pressure
		Q.z is my previous pressure
		rr is the distance between me and my neighbor
		r*(P-Q.z)/rr/N is this interaction's contribution
			to the my velocity. It is the gradient of
			the pressure between me and my neighbor
		length(V0-U0)-length(V1-U1) is the space contraction between
			the next and the previous state. Because pressure
			is energy per volume. If the volume contracts
			there must be a higher energy density. If
			the space expands, the energy is disipated
			over a larger area.
		I add P because I completely trade pressures 
			with my neighbors. This is why the best
			sampling patterns are like checkerboards.
			every other point in space should interact
			so that there is a feedback system.
			The physical explanation of this is that
			the fluid is kind of like many billiards
			balls. When two billiards balls collide,
			they completely swap energies.
		I divide the outputs by N to average each interaction.
    */
}

vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
   R = RENDERSIZE.xy;
 	vec2 U0 = U - C(U).xy, // get previous state
         U1 = U + C(U).xy; // get next state
 	float P = 0.; Q = C(U0);
 	N = 4.;// checkerboard sampling pattern
    P += Xd (U0,U,U1,Q, vec2( 1, 0) );
 	P += Xd (U0,U,U1,Q, vec2( 0,-1) );
 	P += Xd (U0,U,U1,Q, vec2(-1, 0) );
 	P += Xd (U0,U,U1,Q, vec2( 0, 1) );
 	Q.z = P;
 	
 	// Init and Walls
 	if (int(FRAMECOUNT) <= 1 || Reset != 0.) Q = vec4(0);
    if (U.x < 1.||U.y < 1.||R.x-U.x < 1.||R.y-U.y < 1.) Q.xy *= 0.;
    //This Jet setup will cause pressure to always increase in the system
 	if (length(U-vec2(0,0.5*R.y)) < 4.) {Q.xy= Q.xy*.9+.1*vec2(.8,0);; Q.w = 1.;}
    if (length(U-vec2(R.x,0.5*R.y)) < 4.) {Q.xy= Q.xy*.9+.1*vec2(-.8,0);; Q.w = 1.;}
	return Q; 
 } 


/*  Hello! Let me try to explain how this works.
	It took me a long time to undertand this stuff
	so bare with me.
	
	Lets start with the ideal gass law:
	PV = nRT
	Temperature is actually the kinetic energy of particles
	So without really getting into it, nRT is the total KE of the system
	So that means PV equals the kinetic energy of a system

	We could arrive at this conclusion in another way as well

	P = F/A         (Pressure is force per area)
	P = F/A * d/d   (here we multiply by 1 in the form of another spacial dimension)
	P = F*d / V     (Area * Length = Volume)
	P = E / V		(Force * distance = Work)
	PV = E  --> P = E/V
	
	So, it follows that pressure is kinetic energy per volume
	
	So lets think of any fluid as a large sum of particles
	that are continuously banging into each other.
	How can you possibly keep track of all that?!
	Well, we can represent different forms of kinetic energy.
	Specifically chaotic kinetic energy and ordered kinetic energy

	so we have 4 channels to work with, so we could have :

r(x)	x direction KE
g(y)	y direction KE
b(z)	and chaotic KE (which is pressure)
a(w)	ink to move around

	Since the medium is moving though itself, any time we
	want to find out about one of these values, we
	have to look them up where they used to be.

	So if I am a pixel and I want to know about my state
	I have to look up the velocity where I am and then
	subtract it from my current location to make a good
	guess of where I was last. In code this looks like this:

	U = my position
	T(U).xy = velocity where I am now
	U - T(U).xy = where I probably was last
	
	Now we want to find out how my pressure (which is
	also my internal kinetic energy) has changed since
	the last state of the fluid.

	Since the last time we knew the state of the
	fluid, each part of the fluid as interacted with its
	neighbors. When one region of space interacts with
	another, it is kind of like they collided with each
	other. Their energies talk to each other in three
	significant ways : 
	
1.	Ordered energy is lost to chaos: A change in 
		the separation between two regions
		corresponds to a linear change in volume.
		A linear change in volume corresponds to a
		linear change in pressure.

2.  The two regions completely exchange pressure.
	 	Think of this as two billiards balls smacking 
		together and trading energies. Or like energy
		traversing a newton's cradle. Remember, pressure
		is actually kinetic energy per volume.

3.  Chaotic energy is converted into ordered energy.
		The space accelerates in the direction of the
		gradient of pressure. Theres are a lot of ways
		to think about this : 
	
			a. A particle sliding down an energy wave

			b. The gradient of energy is a force

			c. Pressure is like pushing outwards, 
				if you have an uneven push, there 
				will be an acceleration. 

	
Conclusion : 
	
-	The space moves each frame according to the ordered
	energy in the space. Ordered energy makes change!

-	Some ordered energy is lost to chaos when the the
	geometry of the space is strained.

-	Directional changes in chaos turn into ordered energy

-	Chaos is continuously traded around like a hot potato
	
	
See Buffer A,B,C, or D and see the function "X"
	to see how this looks in code!



:D Wyatt

*/


vec4 renderMainImage() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;

    U = U/RENDERSIZE.xy;
    vec4 g = texture(BuffA,U,1.);
   	C.xyz = cos(.5-3.*g.w*vec3(1,2,3));
    C = C*sqrt(max(C,0.));
	return C; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderPassC();
	}
	if(PASSINDEX == 3){
		return renderPassD();
	}
	if(PASSINDEX == 4){
		return renderMainImage();
	}
}