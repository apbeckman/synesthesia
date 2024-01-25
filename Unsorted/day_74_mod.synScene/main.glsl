//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********


#define TIME ((TIME*0.3) - 7.)

#define SPEED 0.3
#define FL_H 0.3


#define mx (10.)
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
#define pal(a,b,c,d,e) ((a) + (b)*sin(6.28*((c)*(d) + (e))))
vec3 glowB = vec3(0);

vec3 reflAtten = vec3(1);


vec3 path (float z){
    z *= 0.25;
	return vec3(
    	sin(z + (_uvc.y*FOV+1.0*(cos(z*0.7)+1.0)))*0.7,
    	cos(z + (_uvc.x*FOV+1.0*(cos(z*1.2))+1.0))*(0.6),
        0.
    )*2.;
}

#define pmod(p,x) mod(p,x) - x*0.5
float map(vec3 p){
	float d = 10e6;
    
    // w is used for the lines 
    // and p is used for the tunnel
    
    vec3 w = p; 
    w = abs(w);
    
    // the tunnel is made by the next two lines, otherwise it's just to planes
	p -= path(p.z);
    vec2 try = vec2 (cos(smoothTimeC), sin(smoothTimeC));
    p.xy *= rot(
        sin(w.z*2.9 + p.z*0.7 + sin( w.x*2. + w.z*4. + smoothTimeC*0.15 + 0.5) + w.z*0.1)*(1.6)
    ); 
    //+pow(basshits*0.122*Twitch*syn_Intensity, 2.0)
    float flTop =(-p.y + FL_H )*0.13;
    float flBot =(p.y + FL_H )*0.3;
    float floors = min(flBot, flTop);
    d = min(d,floors);
    d+= WallHeight*0.125;
    float sep = 0.2*(1.0+pow(sin(smoothTimeB*0.1), 2.)*0.1+0.1+separation); // seperation between glowy lines
    
    w.y = pmod(w.y,(sep));
    
    float flash = pow(syn_HighLevel*0.35+syn_MidHighLevel*0.35+syn_Hits*0.25+syn_HighHits*0.125, 2.);
    vec3 z = p;
    // random attenuation to feed to the glowy lines
    float atten = pow(abs(sin(z.z*0.2 + (smoothTimeB*0.1))), 50.);
    float attenC = pow(abs(sin(z.z*0.1  + sin(z.x + smoothTimeB)*0.1 + + sin(z.y*3.)*4. + smoothTimeB*0.1)), 100.);
    float attenB = pow(abs(sin(w.z*0.2  + sin(w.x + smoothTimeB)*0.1 + sin(w.y*0.7)*4. + (w.y*20.) + (smoothTimeB*0.1))), 10.);
    //vec3 col = pal(0.1,0.6 - attenC*0.65,vec3(1.7  - atten*0.4,1.1,0.8),0.2 - atten*0.34 ,0.5 - attenB*0.56 );
    vec3 col = pal(0.16,0.16 - attenC*0.65,vec3(1.7  - atten*01.4,1.61,01.8),0.2 - atten*0.34 ,0.5 - attenB*0.56 );
    col = max(col, 0.0125*(1.0+flash*0.75));

    
    float sc = 60. - atten*65.;
    
    // distance to the glowy lines
    float dGlowzers = max(floors,-abs(w.y) + sep*0.5) - 0.02;
       /*
        if(Flash == 0.) {
    glowB += exp(-dGlowzers*(70.))*reflAtten*col*40.;
    }

    else {*/
    glowB += exp(-dGlowzers*(70.))*reflAtten*col*40.*(0.5+0.5*flash);
   // }
    

    // glow
   // glowB += exp(-dGlowzers*(70.))*reflAtten*col*40.;
    d *= 0.65+Twitch*basshits;
    return d;
}

float march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	float d = 10e6;
	p = ro; t = 0.; hit = false;
    for (int i = 0; i < 150 ; i++){
    	d = map(p);
        //glow += exp(-d.x*20.)*0.001;
        if(d < 0.001){
        	hit = true;
            break;
        }
    	t += d;
        p = ro + rd*t;        
    }

    return d;
}

vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
	vec3 dir = normalize(lookAt - ro);
    vec3 right = normalize(cross(vec3(0,1,0), dir));
    vec3 up = normalize(cross(dir, right));
	return normalize(dir + right*uv.x + up*uv.y);
}


vec3 getNormal(vec3 p){
	vec2 t = vec2(0.001,0);
    return normalize(
    	vec3(
        	map(p - t.xyy) - map(p + t.xyy),
        	map(p - t.yxy) - map(p + t.yxy),
        	map(p - t.yyx) - map(p + t.yyx)
        )
    );
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    uv *= 1. - dot(uv,uv)*-0.2;
    
   // uv.xy *=1.0+(_uvc*PI*FOV*uv.xy)*0.5;
    //uv.xy *= rot(0.1)
    vec3 col = vec3(0);
    
    vec3 ro = vec3(0);
    
    ro.z += mx*2.;
    ro.z += bass_time;

    ro += path(ro.z);
    
    vec3 lookAt = vec3(0,0,ro.z + 1.);

    lookAt += path(lookAt.z);
    //lookAt.xy += (_uvc*PI*FOV*lookAt.xy);
   // uv *=1.0+ WarpFurther*(-1.+PI*uv*PI);
    vec2 xy_noise = vec2(_noise(vec2(sin(TIME*0.2), cos(TIME*0.2))));
    vec3 rd = getRd(ro, lookAt, (uv)+_uvc*PI*FOV);
    rd.xz = _rotate(rd.xz, LookXY.x*PI+_uvc.x*PI*WarpFurther.x+ xy_noise.x*0.2);
    rd.xy = _rotate(rd.xy, Rotate*PI);
    rd.yz = _rotate(rd.yz, LookXY.y*PI+_uvc.y*PI*WarpFurther.y*0.75 + xy_noise.x*0.2);
    //rd.xy *= rot(sin(smoothTime)*0.05);
    
    bool hit; float t; vec3 p;
    
    float bounce;
    
    float firstT = 0.;
    float d;
    for(int i = 0; i < 2 ; i++){
        d = march(ro, rd, p, t, hit);
        vec3 n = getNormal(p);
        if(i == 0)
        
        reflAtten *= vec3(0.48);
        
        rd = reflect(rd, n);
        ro = p + rd*0.2;
    }
    
    
    glowB = max(glowB, 0.);
    
    col += glowB*0.0006;
    
    col = clamp(col,0., 1.);
    col = pow(col, vec3(1.3));
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 



// I fixed up the shader a bit, compared to the original and added tome color toning

// It is basically is just two perpendicular planes which are rotated depending on the position of the viewer.
// Materials are reflective
// Then some glowy lines are added on top



vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    
    // Radial blur and chromatic abberation
    float steps = 24.;
    float scale = 0.00 + pow(length(uv - 0.5),3.4)*0.1;
    float chromAb = pow(length(uv - 0.5),1.)*2.5;
    vec2 offs = vec2(0);
    vec4 radial = vec4(0);
    for(float i = 0.; i < steps; i++){
    
        scale *= 0.99;
        vec2 target = uv + offs;
        offs -= normalize(uvn)*scale/steps;
    	radial.r += texture(BuffA, target + chromAb*1./RENDERSIZE.xy).x;
    	radial.g += texture(BuffA, target).y;
    	radial.b += texture(BuffA, target - chromAb*1./RENDERSIZE.xy).z;
    }
    radial /= steps;
    
    fragColor = radial*1.5; 
    
    // mimap glow
    //fragColor += texture(BuffA,uv, 6.)*0.1;
    
    // color correction
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.14); 
    
    fragColor *= 1.;
    fragColor.g *= 1.1;
    fragColor.r *= 0.95 + uvn.x*0.7;
    fragColor.g *= 0.95 + uvn.y*0.3;
    fragColor = max(fragColor, 0.);
    // vignette
    fragColor = pow(fragColor, vec4(0.545 + dot(uvn,uvn)*2.)); 
	
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