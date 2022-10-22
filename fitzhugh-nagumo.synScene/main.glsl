vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********

/*
	A Fitzhugh-Nagumo reaction-diffusion system. 
	See this paper for additional information: 
		
	http://arxiv.org/pdf/patt-sol/9401002.pdf

	A large timestep is used to make the system evolve at an interactive rate when limited to 60 FPS.
    The system is unstable using a large timestep with simple Euler integration, so instead it is 
    updated with an exponentially-weighted moving average of the gradient (with time constant tc).
*/

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    const float _K0 = -20.0/6.0; // center weight
    const float _K1 = 4.0/6.0; // edge-neighbors
    const float _K2 = 1.0/6.0; // vertex-neighbors
    const float timestep = 0.7;
    const float a0 = -0.1;
    const float a1 = 2.0;
    const float epsilon = 0.05;
    const float delta = 4.0;
    const float k1 = 1.0;
    const float k2 = 0.0;
    const float k3 = 1.0;
    const float tc = 0.8;

    vec2 mouse = iMouse.xy / RENDERSIZE.xy;
    vec2 vUv = fragCoord.xy / RENDERSIZE.xy;
    vec2 texel = 1. / RENDERSIZE.xy;
    
    // 3x3 neighborhood coordinates
    float step_x = texel.x;
    float step_y = texel.y;
    vec2 n  = vec2(0.0, step_y);
    vec2 ne = vec2(step_x, step_y);
    vec2 e  = vec2(step_x, 0.0);
    vec2 se = vec2(step_x, -step_y);
    vec2 s  = vec2(0.0, -step_y);
    vec2 sw = vec2(-step_x, -step_y);
    vec2 w  = vec2(-step_x, 0.0);
    vec2 nw = vec2(-step_x, step_y);

    vec4 uv =    texture(BuffA, vUv);
    vec4 uv_n =  texture(BuffA, vUv+n);
    vec4 uv_e =  texture(BuffA, vUv+e);
    vec4 uv_s =  texture(BuffA, vUv+s);
    vec4 uv_w =  texture(BuffA, vUv+w);
    vec4 uv_nw = texture(BuffA, vUv+nw);
    vec4 uv_sw = texture(BuffA, vUv+sw);
    vec4 uv_ne = texture(BuffA, vUv+ne);
    vec4 uv_se = texture(BuffA, vUv+se);

    // laplacian of all components
    vec4 lapl  = _K0*uv + _K1*(uv_n + uv_e + uv_w + uv_s) + _K2*(uv_nw + uv_sw + uv_ne + uv_se);

    float a = uv.x;
    float b = uv.y;
    float c = uv.z;
    float d = uv.w;
    
    float d_a = k1*a - k2*a*a - a*a*a - b + lapl.x;
    float d_b = epsilon*(k3*a - a1*b - a0) + delta*lapl.y;
	c = tc * c + (1.0 - tc) * d_a;
	d = tc * d + (1.0 - tc) * d_b;

    a = a + timestep * c;
    b = b + timestep * d;
    
    if (iMouse.z > 0.0) {
    	float mLen = length(iMouse.xy - fragCoord.xy);
    	a += exp(-mLen * mLen / 100.0);
    }
    
    // initialize with noise
    if(FRAMECOUNT<10) {
        fragColor = -0.5 + texture(image16, vUv);
    } else {
        fragColor = clamp(vec4(a, b, c, d), -1., 1.);
    }
    

	return fragColor; 
 } 


/*
	The reaction-diffusion system is visualized with a slightly modified version of 
    Shane's Bumped Sinusoidal Warp shadertoy here:

	https://www.shadertoy.com/view/4l2XWK
    
	The x channel of Buffer A, containing the reaction-diffusion system components,
    is used for the bump mapping function.
*/


// Bump mapping function.
float bumpFunc(vec2 p){ 
    return 0.5 * (texture(BuffA, p).x + 1.0);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    // Screen coordinates.
	//vec2 uv = (fragCoord - RENDERSIZE.xy*.5)/RENDERSIZE.y;
    vec2 uv = fragCoord.xy/RENDERSIZE.xy;
    
    // VECTOR SETUP - surface postion, ray origin, unit direction vector, and light postion.
    vec3 sp = vec3(uv, 0); // Surface posion. Hit point, if you prefer. Essentially, a screen at the origin.
    vec3 rd = normalize(vec3(uv - 1.0, 1.)); // Unit direction vector. From the origin to the screen plane.
    vec3 lp = vec3(cos(TIME/2.0)*0.5, sin(TIME/2.0)*0.5, -1.); // Light position - Back from the screen.
	vec3 sn = vec3(0., 0., -1); // Plane normal. Z pointing toward the viewer.

    vec2 eps = 2.0 / RENDERSIZE.xy;
    
    float f = bumpFunc(sp.xy); // Sample value multiplied by the amplitude.
    float fx = bumpFunc(sp.xy-vec2(eps.x, 0.0)); // Same for the nearby sample in the X-direction.
    float fy = bumpFunc(sp.xy-vec2(0.0, eps.y)); // Same for the nearby sample in the Y-direction.
   
 	// Controls how much the bump is accentuated.
	const float bumpFactor = 0.02;
    
    // Using the above to determine the dx and dy function gradients.
    fx = (fx-f)/eps.x; // Change in X
    fy = (fy-f)/eps.y; // Change in Y.
    sn = normalize( sn + vec3(fx, fy, 0)*bumpFactor );           
   
    
    // LIGHTING
    //
	// Determine the light direction vector, calculate its distance, then normalize it.
	vec3 ld = lp - sp;
	float lDist = max(length(ld), 0.001);
	ld /= lDist;

    // Light attenuation.    
    float atten = min(1./(0.25 + lDist*0.5 + lDist*lDist*0.05), 1.);
    
    atten *= f*f*.5 + .5;

	// Diffuse value.
	float diff = max(dot(sn, ld), 0.);  
    // Enhancing the diffuse value a bit. Made up.
    diff = pow(diff, 2.)*0.66 + pow(diff, 4.)*0.34; 
    // Specular highlighting.
    float spec = pow(max(dot( reflect(-ld, sn), -rd), 0.), 8.); 

    vec3 texCol = texture(image3, sp.xy).xyz;
    
    // FINAL COLOR
    // Using the values above to produce the final color.   
    vec3 col = (texCol * (diff*vec3(1, .97, .92)*1.3 + 0.5) + vec3(1., 0.6, .2)*spec*1.3)*atten;

    // Done. 
	fragColor = vec4(min(col, 1.), 1.);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}