//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 



#define T (script_time*0.15)+(TIME*0.1)+(smooth_basstime*0.1)

// wrighter here 
// -- plant -- //
// fractal in the middle
// mist 
// floating object inside walls 
// floating mysterius tentacles

// everything is going to be glowy, transparent
// except of the floating objects inside of the walls 


// if i have time
// transparent
// polar repetetion, 4 times

vec3 path (float z){
    z *= 0.5;
    return vec3(sin(z), cos(z),0);
}

vec3 glow = vec3(0);

#define pal(a,b,c,d,e) ((a)+(b)*sin((c)*(d) + (e)))

#define rot(j) mat2(cos(j),-sin(j),sin(j),cos(j))

#define pi acos(-1.)
#define tau (2.*pi)

#define kick (floor(T) + pow(fract(T),7.))


vec2 map(vec3 p){
		vec2 d = vec2(10e6);
        float dTrans = 10e7;

        p.y -= 1.;

    	p -= sin(p)/1.;
    
        //p.z += T; 
        vec4 q = vec4(p,1.);
        for(int i = 0; i < 5; i++){

            float dpp = dot(q.xyz,q.xyz);
            q /= dpp;

            q.xyz = abs(q.xyz) - vec3(0.03,0.02,0.4);
			
            //q.xz *= rot(0.5);
            q.xy *= rot(0. + sin( (T + kick*pi/2.)*0.6)*0.14 );
            
            //q.zx *= rot(-0.2 );

        }

        dTrans = min(dTrans, length(q.xyz)/q.w);

        dTrans = min(dTrans, length(q.xy)/q.w);

        dTrans -= 0.003+ Thiccness;

        dTrans *= 0.5; 
        dTrans = abs(dTrans) + 0.004;


        //dTrans = min(dTrans, length(p)-0.4);

        glow += 0.7/(0.01 + dTrans*dTrans*2000.)*pal(0.5,0.5,vec3(0.7,0.2,0.1),1.,length(p) + T );
		
        d.x = min(d.x, dTrans);

        dTrans = 10e8;
        dTrans = min(dTrans, length(q.xz)/q.w);
    	glow += 8.7/(0.01 + dTrans*dTrans*200.)*pal(0.5,0.5,vec3(0.89,0.2,0.1),1.,length(p) + smooth_hightime )*pow(abs(sin(length(p)*4. + T - kick*4.)),40.);

        dTrans = abs(dTrans*0.5) + 0.004;

        d.x = min(d.x, dTrans);

        q.yz *= rot(.0);

    	q.x +=.0;
    	q.z -=.0;
    	
    	q.z -= 0.;
    	q.xy *= rot((0.3)*pi);
    	q = abs(q);
        float db = max(q.z,q.x)/q.w*0.4;
        db -= 0.04;
        //glow += 0.7/(0.01 + dTrans*dTrans*200.)*pal(0.5,0.5,vec3(0.89,0.2,0.1),1.,length(p) + T )*pow(abs(sin(p.x*4. + T)),40.);

        d.x = min(d.x,db);
		
    	float dbb = abs(length(p)-1.5)+0.003; 
        d.x = min(d.x,dbb);
    	glow += exp(-dbb)*0.1;
		
    
    	d.x *= 0.9;
	return d;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    vec3 col = vec3(0.04,0.0,0.0);
    //col = vec3(0.7,0.1,0.5);

    uv *= 1. + dot(uv,uv);
    
    vec3 ro = vec3(0,0,0.0);
    ro.y += sin(T);
    ro.xz += vec2(cos(T)*0.4 + 1.,sin(T)*0.5-1.2)*9.;

    vec3 lookAt = vec3(0);

    lookAt.z = ro.z + 2.;

    lookAt -= lookAt;

    ro += path(ro.z);

    lookAt += path(lookAt.z);


    vec3 dir = normalize(lookAt - ro);

    vec3 right = normalize(cross(vec3(0,1,0), dir));
    vec3 up= normalize(cross(dir, right));

    vec3 rd = normalize(dir + right*uv.x + up*uv.y);


    vec3 p = ro; float t = 0.; vec2 d = vec2(10e7);
    bool didHit = false;



    p.z -= sin((smooth_basstime*0.4) +(TIME*0.3)+ 4.);

    for(int i = 0; i < 90; i++){
        d = map(p);
        if(d.x < 0.0025){
        	didHit = true;
            break;
        }
        
        t += d.x;
        p = ro + rd*t;
    }

    col += glow*0.04;

    
    
    
    if(didHit){
    	vec2 t = vec2(0.0025,0);
        vec3 n = normalize(map(p).x - vec3(map(p-t.xyy).x,map(p-t.yxy).x,map(p-t.yyx).x));
        
        
        #define AO(j) clamp(map(p + n*j).x/j,0.,1.)
        
        col += 2.*AO(0.9)*AO(0.1);

    }

    col = 1. - col;
    col = smoothstep(0.,1.,col);
    
    col *= 1. - dot(uv,uv)*0.1;
    
    //col = mix(col,vec3(0.5),smoothstep(0.,1.,t*0.001));
    //col 
    /*
    col /= 1. - col*0.2;
    
    */
    //col = smoothstep(0.,1.,col*1.);
    
    col = pow(col,vec3(0.454545));

    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}