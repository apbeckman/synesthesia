

			//******** BuffA Code Begins ********


#define dmin(a,b) a.x < b.x ? a : b
#define pi acos(-1.)
#define tau (2.*pi)

#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))

#define pal(a,b,c,d,e) ((a) + (b)*sin((c)*(d) + (e)))

#define pmod(p,z) mod(p,z) - 0.5*z

vec3 glow = vec3(0);

float sdBox(vec3 p, vec3 s){
	p = abs(p) - s;
    return max(p.x,max(p.y,p.z));
}	


vec3 qG;
vec2 idG;
vec2 map(vec3 p){
	vec2 d = vec2(10e8);

    vec3 q = p;
    
    float modD = 1.;
    float reps = 10.;
    
    q = vec3(reps*atan(q.z,q.x)/tau ,log(length(q.xz)) ,q.y);
    
    //q.
    
    //q.z += length(p.xz)*0.6;
    //q.z += exp(length(p.xz)*0.4)*0.1;
    float h = exp(length(p.xz)*0.7)*0.2;
    q.z -= h;
    q.z += exp(-length(p.xz)*4.)*0.7;
    //q.y *=  1. - exp(-length(p.xz)*5.)*6.4;
    q.y *= reps/tau;
    q.y += TIME;
    
    float id = floor(q.y);
    //q *= 1.;
    
    q.x += TIME*0.4 + id*0.6 + q.y*0.1;
    
    //q.xy *= rot(0.1 + id);
    idG = floor(q.xy/modD);
    q.xy = pmod(q.xy, modD);
    
    qG = q;
    
    //q.xz = sin(q.xz);
    float bW = 0.6;
    
    //float dB = length(q) - 0.3;
    //for()
    
    q.xy *= rot(0.25*pi);
    
    float dB = sdBox(q, vec3(bW*2.,bW*2.,bW));
    //dB = max(dB, -sdBox(q, vec3(0.1, bW*0.2,0.7)));
    float s = 0.2*bW;
    dB = min(dB, sdBox(q, vec3(bW*s, bW*s,bW*1.2)));
    
    s *= 2.;
    
    q.xy *= rot(0.5);
    float dBb = sdBox(q, vec3(bW*s, bW*s,bW*1.1));
    d = dmin(d,vec2(dBb, 1.) );
    
    
    
    d = dmin(d, vec2(dB, 2.));
    
    
    s = 0.2*bW;
    q.xy *= rot(0.5*pi);
    float dC = sdBox(q, vec3(bW*s*1., bW*s*1.,bW*1.3));
    d = dmin(d, vec2(dC, 3.));
    
    
    
    q = abs(q);
    //q.xz -= 0.;
    //q.yz *= rot(0.25*pi);
    
    q.z -= 0.4;
    float dD = sdBox(q, vec3(bW*s*1., bW*s*1.,bW*20.5));
    
    
    //glow += 0.8/(0.01 + dD*dD*20.)*pal(0., 0.9,vec3(0.9,0.4,0.8), 0.7,1.2);
    p.y -= h;
    glow += 1.8/(0.06 + dD*dD*200.)*sin(vec3(0.8,0.5,0.1) + vec3(0,-0. - length(p.xz)*0.1,0))*pow(abs(sin(TIME + idG.y)), 20.)* smoothstep(1.,0.,length(p.y)*0.6)* pow(smoothstep(0.,1.,length(p.xz)*4.2),7.); // pow(smoothstep(0.,1.,length(p.xz)*1.2),7.);
    
    dD = abs(dD) + 0.006 + smoothstep(0.,1.,length(p.y)*0.05);
    d = dmin(d, vec2(dD, 3.));
    //= dmin(d, vec2(length(p)- 0.1, 2.));
    d.x *= 0.4;
    
    d.x *= exp(log(length(p.xz) + 0.1)/1.5);
    
    return d;
}

float dith;

vec2 march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	vec2 d= vec2(10e8);

    t = 0.; hit = false; p = ro;
    
    for(int i = 0; i < 150; i++){
    	d = map(p);
        d.x *= dith;
        
        if(d.x < 0.001){
        	hit = true;
            break;
        }
        t += d.x;
        p = ro + rd*t;
    }
    
    return d;
}

vec3 getNormal(vec3 p){
	vec2 t = vec2(0.0001,0.);
    return normalize(
    	vec3(
        	map(p + t.xyy).x - map(p - t.xyy).x,
        	map(p + t.yxy).x - map(p - t.yxy).x,
        	map(p + t.yyx).x - map(p - t.yyx).x
        )
    );

}

vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
	vec3 dir = normalize(lookAt - ro);
	vec3 right = normalize(cross(vec3(0,1,0), dir));
	vec3 up = normalize(cross(dir, right));
    return normalize(dir + right*uv.x + up*uv.y);
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    uv *= 1. + dot(uv,uv)*0.1;
    
    vec3 col = vec3(0);
    dith = mix(0.96,1.,texture(image30, RENDERSIZE.x*(uv/256.)).x);
    
    vec3 ro = vec3(2.501,3,0.001);
    vec3 lookAt = vec3(0.05);
    
    lookAt.xz += vec2(sin(TIME)*1.2,cos(TIME)*0.5)*0.4;
	vec3 rd = getRd(ro, lookAt, uv); 
    vec3 p;float t; bool hit;
    
    vec2 d = march(ro, rd, p, t, hit);
    
    
    if(hit){
    	vec3 n = getNormal(p);
        vec3 l = normalize(vec3(1));
        
        
        #define AO(a) clamp(map(p + n*a).x/a,0.,1.)
        float dnrd = max(dot(n, - rd),0.);
        
        float diff = max(dot(n, l),0.);
        float spec = pow(max(dot(n, normalize(l - rd)),0.), 20.);
        float fres = pow(1.- dnrd, 5.);
        
    	//col += (0.5 + n*0.5)*2.;
        
        vec3 q = qG;
        
        float modD = 0.1;
        q = abs(q);
        
        q.xy *= rot(0.25*pi);
        q.y -= 0.5;
        
        q = abs(q);
        q.xy *= rot(0.25*pi);
        q = abs(q);
        q.x -= 0.5;
        q.xy *= rot(-0.25*pi);
        
        
        float ao = mix(0.4,AO(0.1)*AO(0.2)/0.26,smoothstep(0.,1.,length(p.xz)*1.));
        //float ao = AO(0.2)*AO(0.4)*AO(0.3)/0.16;
        //float ao = AO(0.2)*AO(2.4)*AO(0.8)/0.3*3.;
        float idc = floor(max(q.x,q.y)/ modD + TIME*4.*0.25);
        
        
        
        if(d.y == 2. || d.y == 1.){
            //col += pal(0.5,1., vec3(.7,0.1,-0.3),0.8,idc*1.5 + idG.y*0.5 + TIME*0.1);
            col += pal(0.,0.9, vec3(.7,0.1,-0.3),0.8,idc*1.5 + length(p.xz)*1. + TIME*1.);
		

            col = max(col, 0.);
        } else if(d.y == 3.) {
            //col += pal(0.5,1., vec3(.7,0.2,0.),0.8,idc*1.5 + idG.y*0.5);
            col += 0.4*fres;

            
            
            //col += shiftHue(col, 0.98).xyz;
        }
        //col += pal(0.5,0.5, vec3(.1,0.9,1.8),0.8,idc*1.5 );
        if(d.y == 1.){
        	//col = sin(col*10.);
        }
        
        
        //col *= 0.2+ ao;
        col *= ao*1.4;
        
        col *= smoothstep(0.,1.,length(p.xz)*1.);
        
        //q.xz
        
        //col += qG;
        
    	    
    }
    
    
    col += glow*0.005;
    
    //col = pow(col, vec3(0.454545));

    col = mix(col,vec3(0.06,0.02,0.05), smoothstep(0.,1.,t*0.04));
    
    //col *= 1. - dot(uv,uv)*0.7;
    
    // Output to screen
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


// THIS IS LOGPOLAR MAPPING. THERE IS A REALLY NICE ARTICLE ON HOW TO DO IT. FUN STUFF.
// https://www.osar.fr/notes/logspherical/

// glow from here https://www.shadertoy.com/view/XlBSRz
// pallete - inigo quilez
// dither, ao - nusan
// 

// radiual chromab in this buffer

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    
	fragColor = vec4(0);
    // Radial blur
    float steps = 30.;
    float scale = 0.00 + pow(dot(uvn,uvn),5.5)*0.4;
    float chromAb = pow(length(uv - 0.5),2.9)*0.1;
    vec2 offs = vec2(0);
    vec4 radial = vec4(0);
    for(float i = 0.; i < steps; i++){
        scale *= 0.97;
        vec2 target = uv + offs;
        offs -= normalize(uvn)*scale/steps;
    	radial.r += texture(BuffA, target + chromAb*1.4/RENDERSIZE.xy).x;
    	radial.g += texture(BuffA, target).y;
    	radial.b += texture(BuffA, target - chromAb*1./RENDERSIZE.xy).z;
    }
    radial /= steps;
    
    fragColor += radial;
    
    //fragColor = 1. - fragColor;
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.6);
    //fragColor = pow(fragColor, vec4(1. + dot(uvn,uvn))*1. );
    fragColor = max(fragColor, 0.);
    //fragColor.b *= 1. + uv.x*0.4;
    //fragColor *= 1. - dot(uvn,uvn)*2.;
    fragColor *= 1. - dot(uvn*0.7,uvn*0.7)*1.;
    fragColor = max(fragColor*1.1, 0.);
    fragColor = pow(fragColor, vec4(0.4545 ));
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