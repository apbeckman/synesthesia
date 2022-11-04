//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********


#define dmin(a,b) a.x<b.x ? a : b

// first rot is broken!!
#define rot(x) mat2(cos(x),sin(x),sin(x),cos(x)) 
#define rotgood(x) mat2(cos(x),-sin(x),sin(x),cos(x))
#define pmod(p,x) mod(p,x) - 0.5*x

float smoothTime = (smooth_basstime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.75;
float smoothTimeB = (smooth_hightime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.95;
float smoothTimeC = (smooth_midtime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.95;

vec3 glow = vec3(0);
vec3 glowb = vec3(0);
vec3 att = vec3(1);
float side = 1.;

float sdBox(vec3 p, vec3 s){
	p = abs(p) - s;
    return length(max(p,0.0)) + min(max(p.x, max(p.y,p.z)),0.);
}
#define pi (acos(-1.))


vec2 o(vec3 p){
    vec2 d = vec2(10e7);
    
    //p.xy *= rot(0.25*pi);
    p = abs(p);
    float dBox = sdBox(p, vec3(0.1,0.5,0.7));
    p.x -= .45;
    float dBoxb = sdBox(p, vec3(0.08,0.8,0.5));
    
    p = abs(p);
    
    p.y -= 0.45;
    float dBoxc = sdBox(p, vec3(0.16));
    
    p.y -= 0.4;
    p.z -= 0.7;
    p = abs(p);
    
    p.x -= 0.05;
    p.z -= 0.15;
    float dBoxg = max(p.x,p.z);
	
    d = dmin(d, vec2(dBox, 5.));
    d = dmin(d, vec2(dBoxb, 21.));
    d = dmin(d, vec2(dBoxc, 9.));;
    d = dmin(d, vec2(dBoxg, 21.));
    
    return d;
}
float mmmm;
vec2 map(vec3 p){
	vec2 d = vec2(10e6);

    
    vec3 q = p;
    float modD = 9.;
    
    vec2 id = floor(p.xz/modD);
    p.xz = pmod(p.xz,modD);
    
    
    float mt = sin(((smoothTimeB*0.25)+crystalrate)*1.275 + id.x + sin(id.y)*0.4 + id.y*0.4);
    mmmm = pow(abs(mt), 0.5-mod1)*sign(mt)*(0.2125+crazyFloor);
    
    //glowb += exp(-length(p)*5.)*2.*vec3(0.3,0.9,0.9);
    
    vec3 cc = vec3(0.7 + mmmm*2.,.9,0.9);
    glowb += exp(-length(p)*(5. - mmmm))*2.*cc*att;
    
    
    for(float i = 0.; i< 4.; i++){
    	p = abs(p);
        p.xy *= rot((0.5)*pi);
        p.x -= 0.8+ mmmm;
        p.t -= 1.5 ;
        p.zy *= rot(0.25*pi);
        p.z -= 1.;
        //p.xz *= rot(0.5*pi);
        
    }
    
    //p.x -= sin(TIME);
    
    vec2 dO = o(p);
    
    float dBall = length(p) - 0.1;
    
    d = dmin(d, dO);
    
    q = abs(q);
    q.y -=6.;
    
    d = dmin(d, vec2(abs(q.y) - 0.1, 21.));
    d.x *= 0.6;
    
    glow += exp(-max(d.x, 0.)*10.)*att*2.;
    
    return d;
}

vec2 march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	vec2 d = map(ro);
    
    if(d.x < 0.2)
        ro += rd*0.08;
    
    hit = false; t = 0.; p = ro;
    
    for(float i = 0.; i< 140.; i++){
    	d = map(p);
        d.x *= side;
        if(d.x < 0.00008){
        	hit = true;
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
	return normalize(dir + (right*uv.x + up*uv.y)*0.8);
}
vec3 getNormala(vec3 p){
	vec2 t = vec2(0.0001,0.);
	return normalize(map(p).x - vec3(
    	map(p - t.xyy).x,
    	map(p - t.yxy).x,
    	map(p - t.yyx).x
    ));
}
vec3 getNormal(vec3 p){
	vec2 t = vec2(0.0001,0.);
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
    
    //rro += smooth_basstime;
    
    ro.z += smoothTime*2.4;
    
    vec3 lookAt = vec3(0,0,ro.z + 1.5);
    
    vec3 rd = getRd(ro, lookAt, uv);
    rd.yz = _rotate(rd.yz, lookXY.y*PI);
    rd.xy = _rotate(rd.xy, lookXY.x*PI);
    
    vec3 p; float t; bool hit; float tA;
	float tF;
    side = sign(map(ro).x);
    vec2 d;
    int iters = 2;
    for(int i = 0; i < iters + min(FRAMECOUNT, 0); i++){
        d = march(ro, rd, p, t, hit);
        
        if(i == 0) 
            tF = t;
        tA = max(t, tA);
        if(hit){
			vec3 n = getNormal(p)*side;
			vec3 l = normalize(vec3(1));
            
            float diff = max(dot(n,l), 0.);
            float spec = pow(max(dot(reflect(l, -rd),n),0.), 10.);
            float fres = pow(1. - max(dot(n,-rd), 0.), 1.);
            #define ao(j) clamp(map(p + n*j).x/j, 0., 1.)
            #define sss(j) smoothstep(0., 1.,map(p + l*j).x/j)
            float a = ao(0.34)*ao(0.1)*ao(0.1);
    				        
            vec3 lCol = vec3(0.3,0.7,1.);
            if(d.y > 10.){
                col += 0.04*(pow(fres, 1.)*lCol + glow*0.0002 + lCol*spec*0.8)*att*a;
                ro = p;
                //rd = refract(rd, n, 0.99);
                //side *= -1.;
                rd = reflect(rd, n);
                att *= vec3(0.6,0.5,0.6)*0.4;
            } else if(d.y > 5.){
                //col += vec3(0.7,0.7,0.4)*att*(a  + glow*0.0004 + sss(0.2 ));
                col += vec3(0.7,0.7,0.4)*att*fres*2.*(glow*0.0004 );
                //rd = reflect(rd, n);
                ro = p;
                //rd = refract(rd, n, 0.95);
                //side *= -1.;
                //att *= 0.8;
            	break;
            } else {
            	col += vec3(0.4,0.2,0.1)*((0.10 + diff*1.*fres)*att + spec*0.4*lCol)*a*att;
            	break;
            }
            
            //#define FOG vec3(0.25,0.14,0.32)*0.06
            //#define FOG vec3(0.15,0.14,0.28)*0.06
            #define FOG vec3(0.15 + mmmm*0.2,0.10,0.28)*0.04
            if (i == iters - 1){
    			col = mix(col, FOG*att, smoothstep(0.,1.,tA*0.015));
            
            }

        }
    
    }
    
    
    col += glowb*0.002;
    
    col = mix(col, FOG*0.06, pow(smoothstep(0.,1.,tF*0.015), 1.4));
    //col = mix(col, FOG*0.06, pow(smoothstep(1.,0.,exp(-(abs(p.y) - 5.)*20.)), 1.4));
    
    //col = mix(col, vec3(0.1,0.54,0.512)*0.06, smoothstep(0.,1.,tA*0.03));
    
    
    

    // Output to screen
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


// thx to evvvvil and nusan for some of the techniques here
// radial blur and chromatic abberation in this buffer
// The whole thing is built around a incorrectly written rot() function xd


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    

    // Radial blur
    float steps = 1.;
    float scale = 0.00 + pow(length(uv - 0.25),1.)*0.1;
    //float chromAb = smoothstep(0.,1.,pow(length(uv - 0.5), 0.3))*1.1;
    float chromAb = pow(length(uv - 0.5),1.)*0.8;
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
    
    fragColor = radial*25.; 
    //fragColor.r *= 1. + uv.x*0.8;
    //fragColor.g *= 1. + uv.y*0.7;
    //fragColor.b *= 1. + uv.y*0.7;
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.4);
    //1fragColor *= 18.;
    
    
    fragColor = max(fragColor, 0.);
    fragColor = pow(fragColor, vec4(0.4545 + dot(uvn,uvn)*1.));
    //fragColor *= 1. - dot(uvn,uvn)*0.6;
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