

// Author: ocb
// Title: Boxlight 3D

// Mixing "Phot'on stage" (https://www.shadertoy.com/view/4sVcz3)
// with "HOPE" (https://www.shadertoy.com/view/MllfDX)
// for a 3D version of "Phot'on stage".
//
// This shader need to be anti-aliased, using the same process as in
// https://www.shadertoy.com/view/ldyyDh (one ray AA)
//
// MOUSE ENABLE
//

// comment this line to avoid anti-aliasing process
#define AA

// Increase this if perfos permit (max voxel distance)
#define MAXSTEP 60

#define PI 3.141592653589793
#define INFINI 1000000.


float H(in float v) { 						
    return fract(sin(v) * 437585.);
}


float H2(in vec2 st) { 						
    return fract(sin(dot(st,vec2(12.9898,8.233))) * 43758.5453123);
}


vec4 getNextPlan(in vec2 xz, in vec2 v){
    vec2 s = sign(v);
    vec2 d = step(0.,s);
	vec2 dtp = (d-fract(xz*.13))/.13/v;
    vec2 dtq = (d-fract((xz-.7)*.17))/.17/v;
    vec2 dtr = (d-fract((xz-.3)*.33))/.33/v;

    vec2 dmin = min(min(dtp,dtq),dtr);
    float tmin = min(dmin.x, dmin.y);
    
    s *= -step(dmin,vec2(tmin));
    
    return vec4(s.x,0.,s.y,tmin);
}


vec3 grid(in vec2 xz, in float t)
{
    t *= .01;
    vec2 p = floor(xz*.13)/.13+t;
    vec2 q = floor((xz-.7)*.17)/.17+.7*t;
    vec2 r = floor((xz-.3)*.33)/.33+1.3*t;
    
    vec3 tex =  .8*texture(image16,p).rgb 
    		  + .6*texture(image16,q).rgb
    		  + .4*texture(image16,r).rgb;
    
    tex *= smoothstep(0.,1.8,tex.r);
    return tex;
    
}

float map(in vec3 tex){ return 6.*tex.r*tex.r; }

float Mmap(in vec3 tex) {return 3.*(tex.r*tex.r		// Music map
    							  *(texelFetch(iChannel1,ivec2(1,0),0).x
                                  + texelFetch(iChannel1,ivec2(400,0),0).x))
    					;}

vec3 setCol(in vec3 tex, in vec2 xz)	// box color
{
    vec2 bp = 
        	 smoothstep(vec2(.4),vec2(.5),abs(fract(xz*.13)-.5))
    		+ smoothstep(vec2(.4),vec2(.5),abs(fract((xz-.7)*.17)-.5))
    		+ smoothstep(vec2(.4),vec2(.5),abs(fract((xz-.3)*.33)-.5));
            
    vec3 color = tex.rgb*tex.rgb;//*texture(image16,xz*.0025).rgb;
    color *= .3*(bp.x+bp.y)+.3;
    color = clamp(color,0.,1.);

    return color;
}

vec3 setRayCol(in vec3 tex, in vec2 xz)	// light ray color
{
    vec2 bp = 3.-smoothstep(vec2(.0),vec2(.5),abs(fract(xz*.13)-.5))
    		- smoothstep(vec2(.0),vec2(.5),abs(fract((xz-.7)*.17)-.5))
    		- smoothstep(vec2(.0),vec2(.5),abs(fract((xz-.3)*.33)-.5));
            
    vec3 color = tex.rgb*tex.rgb;//*texture(image16,xz*.005).rgb;
    color *= (bp.x+bp.y);

    return color;
}


vec4 trace(in vec3 pos, in vec3 ray, inout vec3 litflux, inout vec3 tex, inout float dAA)
{
    float dh = 0.;
    float t = 0.;
    
    dAA = 1.;				// distance to the border of the voxel
    						// use for AA as transparency factor
    
    vec4 wall = vec4(0.);	// wall.xyz is the normal of the wall. wall.w is the t parameter.
    						// wall is the border of the next voxel
    
    for(int i = 0;i<MAXSTEP;i++){	// entering the voxel run
        
        vec3 p = pos+t*ray;
        
        tex = grid(p.xz,TIME);		// get texture at current time.
        float dh = p.y - Mmap(tex);	// height of the box depend on the red channel of the texture
        
        if( abs(dh)<.25){			// if a box is hit
            
			#ifdef AA
            // finding the closest box border (left/right and up/down) for AA
            dAA = min( getNextPlan(p.xz,ray.xz).w, (.5-2.*abs(dh))+abs(sign(dh)-sign(ray.y)));
            #endif
            
            return vec4(wall.xyz,t);	// you are between top and bottom of a box?
        								// so you hit the previous wall
        }
        
        wall = getNextPlan(p.xz,ray.xz);		// find the next wall
        
        float dhdt = map(grid(p.xz,TIME-.02*p.y+.5+.05*H2(floor(p.xz*20.))));
        				// height of the box back in time depending on p.y
        				// to generate the vertical light flux going up with time
        				// + H2() to create a kind of vertical dithering
        
        litflux += .1*wall.w*wall.w					// to round corner of light flux
            		*setRayCol(tex,p.xz)
            		*smoothstep(6.,20.,dhdt)		// intensity depend on box height reached
            		*min(1.,30./t)*step(0.,dh);		// cutoff with distance and below the box

        float th = -(dh-sign(dh)*.25)/(ray.y+sign(ray.y)*.0001);	// find the box floor/ceiling
        th += step(th,0.)*INFINI;	// eliminate neg value
        
        if(th<wall.w){ // if first hit = floor/ceilling of box
			#ifdef AA
            // distance to the border of the top/bottom surface used for AA
            dAA = (wall.w-th)*(.8*abs(ray.y)+.2);
			#endif
            return vec4(0.,-sign(ray.y),0.,t+th);
        }
        
        t+= wall.w+.0001;			// update global t and do again
    }
    
    return vec4(0.,0.,0.,INFINI);
}


vec3 getCamPos(in vec3 camTarget){
    float 	rau = 50.,
            alpha = _mouse.x/RENDERSIZE.x*4.*PI,
            theta = _mouse.y/RENDERSIZE.y*PI+(PI/2.0001);	
    		
    return rau*vec3(-cos(theta)*sin(alpha),sin(theta),cos(theta)*cos(alpha))+camTarget;
}

vec3 getRay(in vec2 st, in vec3 pos, in vec3 camTarget){
    float 	focal = .5;
    vec3 ww = normalize( camTarget - pos);
    vec3 uu = normalize( cross(ww,vec3(0.0,1.0,0.0)) ) ;
    vec3 vv = cross(uu,ww);
	return normalize( st.x*uu + st.y*vv + focal*ww );
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 st = ( fragCoord.xy - .5*RENDERSIZE.xy ) / RENDERSIZE.y;
    float ti = TIME*.1;
    
    // camera def
    vec3 camTarget = vec3(-50.*sin(2.*ti),25.*_mouse.y/RENDERSIZE.y,-30.*cos(3.*ti));    
    vec3 pos = getCamPos(camTarget);
    pos.y = max(pos.y,Mmap(grid(pos.xz,TIME))+1.); // anti-collision
    
    vec3 ray = getRay(st, pos,camTarget);
    	
    float t = 0.;
    vec3 norm = vec3(0.);
    vec3 color = vec3(0.);
    vec3 litflux = vec3(0.);	// color of the light rays
    vec3 tex;
    float dAA = 1.;		// distance to border of object, used for anti-aliasing
        
    vec4 info = trace(pos, ray, litflux, tex, dAA);
    t = info.w;
    norm = info.xyz;
	vec3 p = pos + t*ray;

    if(t<INFINI){
        vec2 h = vec2(0.,Mmap(tex));
        vec2 side = p.xz*norm.y+(p.xy-h)*abs(norm.z)+(p.zy-h)*abs(norm.x);
        color += setCol(tex,p.xz)*(.7+.005*p.y*p.y*p.y);
        color *= min(1.,20./t);
        color = clamp(color,0.,1.);
        
#ifdef AA
        dAA = min(.2*dAA*RENDERSIZE.x/t , 1.);		// anti-aliasing
        color *= dAA;								// dAA = transparency factor
#endif
        
    }
    color += litflux*(1.-max(ray.y,0.));	// adding photon flux (vanishing in the sky)
    
    
#ifdef AA
    // for AA we need to mix the object color (found above) with the background color
    // here we find the background color
	float T = t;			// to save global distance
    float cAA = 1.-dAA;		// complementary color to add
    
    if(dAA<1.){				// if dAA = 1, full object hit, no need for AA
        
        vec4 wall = getNextPlan(p.xz,ray.xz);	// stepping to the next voxel
        float tt = wall.w+.0001;
        
		// adding light flux to the next voxel
        float dh = p.y - Mmap(tex);				
        float dhdt = map(grid(p.xz,TIME-.03*p.y+.5));
        litflux = .1*wall.w*wall.w*setRayCol(tex,p.xz)*smoothstep(6.,20.,dhdt)*min(1.,30./t)*step(0.,dh);
        
        // now we are positionned at the next voxel => tracing again
        info = trace(p+tt*ray, ray, litflux, tex, dAA);
        t = info.w;
        norm = info.xyz;
		T += t+tt;		// updating global distance
        
        if(t<INFINI){	// a box is hit, setting its color
            p += (t+tt)*ray;
            vec2 h = vec2(0.,Mmap(tex));
            vec2 side = p.xz*norm.y+(p.xy-h)*abs(norm.z)+(p.zy-h)*abs(norm.x);
            vec3 colAA = setCol(tex,p.xz)*(.7+.005*p.y*p.y*p.y);
            colAA *= min(1.,20./T);
            colAA = clamp(colAA,0.,1.);
            
			// adding the missing amount of color
            colAA *= cAA;
            color += colAA;
        }
        color += cAA*litflux;	// adding the missing amount of color light ray
    }
#endif    


	fragColor = vec4(color,1.);
	return fragColor; 
 } 

    
    

vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}