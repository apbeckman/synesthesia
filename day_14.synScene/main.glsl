//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


#define pi acos(-1.)

#define tau (2.*pi)
#define pal(a,b,c,d,e) ((a) + (b)*sin((c)*(d) + (e)))

#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))

vec3 glow = vec3(0);




float map(vec3 p, float t){
	float d = 10e7;
	
    vec4 q = vec4(p, 1.);
    
    float id = floor(p.z);
    
    
    
    for(float i = 0.; i <0. ; i++){
    	q.xyz = abs(q.xyz);
        q.xy *= rot(0.125*pi + i);
    
        q -= 0.1;
    }
    
    
    float dpp = dot(q.xyz,q.xyz);
    
    
    //    q = q/dpp;
    
    vec4 j;
    
    //q /= dot(q.xyz,q.xyz);
    vec3 b = vec3(2.4+geometry.x, 0.6+(geometry.y), 1.+geometryZ);
    b.y*=1.0+ sin(smoothTimeC/10)/10;
    b.x*=1.0+ cos(smoothTimeC/10)/10;
    int Iterations = int(iterations);
    for(int i = 0; i <Iterations ; i++){
    	q.xyz = abs(q.xyz);
        q.xy *= rot(0.25*pi);
        q.xz *= rot(0.3);
    	q.x += 0.1;
    }
    
    for(float i = 0.; i < 9.; i++){
    	q.xyz = abs(mod(q.xyz - 0.5*b, b) )- 0.5*b;
        
        float dpp = dot(q.xyz,q.xyz);
        
        
        dpp = clamp(dpp, 0. ,1.54);
        if(i == 2.)
            j = q;
        
        q = q/dpp;
        if(i == 20.)
        	q.xz *= rot(0.7 );
        
    }
    
    q.xyz *= 1.;
    
    float db = length(j.yx)/q.w - 0.01;
    
        
	float da = length(q.yz)/q.w - 0.02;
    
    float sc = 0.5;
    d = min(db,da);
    //d = da;
    d *= 0.5;
    d += smoothstep(1.,0.,t*.5)*0.7;
    d = abs(d) + 0.00123;
    
    //d += exp(-t*4.)*0.7;
    
    
    vec3 c = vec3(1,1.,1.);
    da *= 0.5;
    db *= 0.5;
    da = abs(da) + 0.003;
    db = abs(db) + 0.003;
    glow += 0.9/(0.01 + da*da*1500.)*c;
    glow -= 0.9/(0.01 + db*db*2000.)*c*vec3(0.,0.8,2.);
    return d;
}
vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
	vec3 dir = normalize(lookAt - ro);
	vec3 right = normalize(cross(vec3(0.,1.,0), dir));
	vec3 up = normalize(cross( dir, right));
	return normalize(dir + right*uv.x + up*uv.y);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    vec3 col = vec3(0.9,0.6,0.4);

    vec3 ro = vec3(0);
    ro.z += smoothTime*0.05;
    
    float T = smoothTime*0.05;
    ro.xy += vec2(cos(T), sin(T))*3.;
    
    vec3 lookAt = vec3(0);
    
    lookAt.z = ro.z + 4.;
    
    vec3 rd = getRd(ro, lookAt, uv);
    rd.yz = _rotate(rd.yz, lookXY.y*PI);
    rd.xy = _rotate(rd.xy, lookXY.x*PI);

    float d;
    vec3 p = ro; float t = 0.; bool hit = false;
    
    for(int i = 0; i < 60; i++){
    	d = map(p, t);
        if(d < 0.001){
        	hit = true;
            //break;
        }
		t += d;
    	p = ro + rd*t;
    }
    
    
    col -= glow*0.001;
    
    col = pow(col, vec3(0.454545))*(1.0+highhits*0.8);
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}