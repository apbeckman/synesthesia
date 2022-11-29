//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********


#define dmin(a,b) a.x<b.x ? a : b

#define pi (acos(-1.))
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x)) 
#define pmod(p,x) mod(p,x) - 0.5*x

#define pal(a,b,c,d,e) (a + b*sin(c*d + e))

vec3 glow = vec3(0);
vec3 glowb = vec3(0);
vec3 att = vec3(1);
float side = 1.;

float sdBox(vec3 p, vec3 s){
	p = abs(p) - s;
    return length(max(p,0.0)) + min(max(p.x, max(p.y,p.z)),0.);
}


vec2 o(vec3 p){
    vec2 d = vec2(10e7);
    
    //p.xy *= rot(0.25*pi);
    p = abs(p);
    float dBox = sdBox(p, vec3(0.1,0.5,0.7));
    p.x -= 1.75;
    float dBoxb = sdBox(p, vec3(0.08,0.8,.5));
    
    p = abs(p);
    p.y -= 1.45;
    
    p.xy *= rot(0.5*pi);
    float dBoxc = sdBox(p, vec3(0.16, 0.15,6.5));
    
    
    p.y -= 0.4;
    p.z -= 0.7;
    p = abs(p);
    
    p.x -= 0.05;
    p.z -= 0.15;
    float dBoxg = max(p.x,p.z);
	
    d = dmin(d, vec2(dBox, 5.));
    d = dmin(d, vec2(dBoxb, 21.));
    d = dmin(d, vec2(dBoxc, 7.));;
    d = dmin(d, vec2(dBoxg, 21.));
    
    return d;
}
float mmmm;
vec2 map(vec3 p){
	vec2 d = vec2(10e6);

    vec3 q = p;
    float wave = pow(abs(sin(q.z*0.2+smoothTimeB*.5)), 200.);
    float waveb = pow(abs(sin( 0.75*(q.z*0.2+script_time))), 400.);
    float modD = 32.;
    
    p.xz -= modD*0.5;
    vec2 id = floor(p.xz/modD);
    p.xz = pmod(p.xz,modD);
    p.z -= modD*0.2;
    
    
    float mt = sin(smoothTimeB*1. + id.x + sin(id.y)*0.4 + id.y*0.4);
    mmmm = pow(abs(mt), 0.125)*sign(mt)*0.2;
    
    vec3 cc = vec3(0.7 + mmmm*2.,.9,0.9);
    for(float i = 0.; i< 5.; i++){

        if(i < 3.)
        	if(p.x < p.z) {p.xz = p.zx;}

    	p = abs(p);        
        p.xy *= rot((0.25)*pi);
        p.x -= 1.5 ;
        p.t -= 3.5  ;
        p.z -= 4.5;
        if(i == 3.){
            vec3 z = p;
            z += 0.4;
            z = abs(z);
            z -= waveb*0.2;
            z = abs(z);
            float dp = max(z.x,z.y) - 0.05;
            d = dmin(d, vec2(abs(dp) + 0.01, 0.));
            dp = max(dp, 0.);
            vec3 ccc = pal(1.,1.,vec3(0.8,2.2,2.5),0.5,0.5 + smoothTimeB + q.z);
            
            //glowb += 0.001/(0.001 + dp*dp*6.)*ccc*att;
            //glowb += 0.155/(0.40 + dp*dp*dp*dp*(10000. - wave*10000.))*ccc*att;
            //glowb += 0.155/(0.40 + dp*dp*dp*dp*(10000. - 5000.))*ccc*att;
            glowb += 0.155/(1.80 + dp*dp*dp*dp*(10000. - 5000.))*ccc*att;
       
            //glowb += 0.001/(0.001 + dp*dp*6.)*pal(0.5,0.5,vec3(0.87,0.4,0.5), 0.5,3.5 + 1.*TIME)*att;
        
        }
        
    }
    
    
    float md = 0.23;
    vec3 h = p;
    h = pmod(h, md);
    float dpg = sdBox(h,vec3(md*0.1));
    //d = dmin(d, vec2(abs(dpg) + 0.2, 0.));
    dpg = max(dpg, 0.);
    //glowb += 0.001/(0.01 + dpg*dpg*20.)*pal(0.5,0.5,vec3(0.87,0.4,0.5), 0.5,3.5 + 0.*TIME)*att*wave;
    glowb += 0.001/(0.01 + dpg*dpg*20.)*pal(0.5,0.5,vec3(0.87,0.4,0.5), 0.5,3.5 + 0.*smooth_hightime)*att*wave;
    glowb *= (0.999+highhits*0.0125);

    
    
    
    //p.x += 2.;
    //p.x -= sin(TIME);
    
    vec2 dO = o(p);
    
    float dBall = length(p) - 01.1;
    
    d = dmin(d, dO);
    
    q = abs(q);
    q.y -=6.;
    
    //d = dmin(d, vec2(abs(q.y) - 0.1, 21.));
    d.x *= 1.;
    d.x *= 0.7;
    glow += exp(-max(d.x, 0.)*10.)*att*2.;
    
    return d;
}

vec2 march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	vec2 d = map(ro);
    
    if(d.x < 0.3)
        ro += rd*0.025;
    
    hit = false; t = 0.; p = ro;
    
    for(float i = 0.; i< 150.; i++){
    	d = map(p);
        d.x *= side;
        if(d.x < 0.0001){
        	hit = true;
        }
        if(t > 50.){
        	break;
        }
    	t += d.x;
        p = ro + rd*t;
    }
    return d;
}

vec3 getRd(vec3 ro,vec3 lookAt,vec2 uv){
	vec3 dir = normalize(lookAt - ro);
    vec3 right = normalize(cross(vec3(0,1,0), dir));
    vec3 up = normalize(cross(dir, right));
	return normalize(dir + (right*uv.x + up*uv.y)*(0.6+FOV));
}
vec3 getRdIso(inout vec3 ro,vec3 lookAt,vec2 uv){
	vec3 dir = normalize(lookAt - ro);
    
    vec3 right = normalize(cross(vec3(0,1,0), dir));
    vec3 up = normalize(cross(dir, right));
    
    ro += (right*uv.x + up*uv.y)*6.;
	return dir;
}
vec3 getNormal(vec3 p){
	vec2 t = vec2(0.0001,0.);
	return normalize(map(p).x - vec3(
    	map(p - t.xyy).x,
    	map(p - t.yxy).x,
    	map(p - t.yyx).x
    ));
}
vec3 getNormala(vec3 p){
	vec2 t = vec2(0.0004,0.);
	return normalize(-vec3(
    	map(p - t.xyy).x - map(p + t.xyy).x ,
    	map(p - t.yxy).x - map(p + t.yxy).x ,
    	map(p - t.yyx).x - map(p + t.yyx).x 
    ));
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
	vec2 uvn = uv;
    uv *= 1. + dot(uv,uv)*0.7;
    
    vec3 col = vec3(0);
    
    vec3 ro = vec3(0);
    
    ro.xy += vec2(cos(smoothTime*0.1), sin(smoothTime*0.125)*0.1)*0.24;
    
    
    float mp = 0.25;
    ro.xy -= mp;
    //ro.x += iMouse.x/RENDERSIZE.y*mp;
    //ro.y += iMouse.y/RENDERSIZE.y*mp;
    
    ro.z -= 10.;
    ro.z += smoothTime*2.;
    
    vec3 lookAt = vec3(0,0,ro.z + 3.3);
    lookAt.xy += FOV*_uvc*PI;
    vec3 rd = getRd(ro, lookAt, uv);
    rd.xy *= rot(-cos(TIME*0.75)*0.1);
    
    
    vec3 p; float t; bool hit; float tA;
	float tF;
    side = sign(map(ro).x);
    vec2 d;
    int iters = 1;
    float fres;
    for(int i = 0; i < iters; i++){
        d = march(ro, rd, p, t, hit);
        
        if(i == 0) 
            tF = t;
        tA = max(t, tA);
        if(hit){
			vec3 n = getNormal(p)*side;
			vec3 l = normalize(vec3(1));
            
            float diff = max(dot(n,l), 0.);
            float spec = pow(max(dot(reflect(l, -rd),n),0.), 20.);
            fres = pow(1. - max(dot(n,-rd), 0.), 5.);
            #define ao(j) clamp(map(p + n*j).x/j, 0., 1.)
            #define sss(j) smoothstep(0., 1.,map(p + l*j).x/j)
            
            
            float a = ao(0.4)*ao(0.1)*ao(0.2);
            //a = max(a, 0.4);
    				        
            vec3 lCol = vec3(0.4,0.7,1.);
            if(d.y > 10.){
                vec3 refl = (spec*0.05*vec3(0.2,0.1,0.4) + pow(fres, 1.)*0.09*pal(1.,1.,vec3(0.8,2.2,2.5),0.5,0.5 + dot(n, -rd) *30. + (script_time*0.125)))*att;
               
                col += refl*0.1;
                att *= vec3(0.4,0.5,0.6)*0.1;
            } else if(d.y > 5.){
                
                vec3 refl = (pow(fres, 1.) + 0.005)*0.09*pal(1.,1.,vec3(0.8,2.2,2.5),0.5,0.5 + dot(n, -rd) *30. + TIME)*att;
                col += refl*a;
                ro = p;
            	break;
            } else {
                vec3 refl = pow(fres, 1.)*2.5*pal(1.,1.,vec3(1.8,2.2,2.5),0.5,0.5 + dot(n, -rd) *0.6 + sin(TIME)*0.1 + 3.4)*att;
                col += 6.*refl*((0.10 + diff*1.*fres)*att + spec*0.4*lCol)*att*a;
            	break;
            }
            #define FOG vec3(0.05 + mmmm*0.,0.10,0.28)*0.04
            if (i == iters - 1){
    			col = mix(col, FOG*att, smoothstep(0.,1.,tA*0.015));
            }

        }
    
    }
    col = mix(col, FOG*0.06, pow(smoothstep(0.,1.,tF*0.025 + pow(fres,7.)*2. - 0.0), 1.4));
    
    
    col += glowb*0.002;
    
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


// Fork of "Day 95" by jeyko. https://shadertoy.com/view/wssyDM
// 2020-03-23 18:49:46

// Fork of "Day 94" by jeyko. https://shadertoy.com/view/tdXcWM
// 2020-03-23 10:08:54

// thx to evvvvil and nusan for some of the techniques here
// radial blur and chromatic abberation in this buffer


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    

    // Radial blur
    float steps = 20.;
    float scale = 0.00 + pow(length(uv - 0.5),3.)*0.1;
    //float chromAb = smoothstep(0.,1.,pow(length(uv - 0.5), 0.3))*1.1;
    float chromAb = pow(length(uv - 0.5),2.)*4.3;
    vec2 offs = vec2(0);
    vec4 radial = vec4(0);
    for(float i = 0.; i < steps; i++){
    
        scale *= 0.97;
        vec2 target = uv + offs;
        offs -= normalize(uvn)*scale/steps;
    	radial.r += texture(BuffA, target + chromAb*1./RENDERSIZE.xy).x;
    	radial.g += texture(BuffA, target).y;
    	radial.b += texture(BuffA, target - chromAb*1./RENDERSIZE.xy).z;
    }
    radial /= steps;
    
    fragColor = radial*85.; 
    //fragColor.r *= 1. + uv.x*0.8;
    //fragColor.g *= 1. + uv.y*0.7;
    //fragColor.b *= 1. + uv.y*0.7;
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.4);
    //1fragColor *= 18.;
    
    
    fragColor = max(fragColor, 0.);
    fragColor = pow(fragColor, vec4(0.4545 + dot(uvn,uvn)*0.2));
    fragColor *= 1. - dot(uvn,uvn)*0.9
        ;
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