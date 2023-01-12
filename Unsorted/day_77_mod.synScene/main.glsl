vec4 iMouse = vec4(RENDERSIZE, 0, 0); 


			//******** Common Code Begins ********

#define SPEED (1.4 )

			//******** BuffA Code Begins ********

vec3 glow = vec3(0);

vec4 noise(vec2 u){
	return texture(image30, (u)/512.);
}
#define pmod(p,x) (mod(p,x) - 0.5*x)
#define dmin(a,b) a.x < b.x ? a : b
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
#define pal(a,b,c,d,e) ((a) + (b)*sin(tau*((c)*(d) + e)))


vec3 path(float z){
    z *= 0.14;
	return vec3(
    	sin(z + cos(z)*0.6),
    	cos(z + sin(z*0.8)*0.5),
    	0.
    )*2.9;
}

// Hex code from BigWings! He has a tutorial on them.
float HexDist(vec2 p) {
	p = abs(p);
    float c = dot(p, normalize(vec2(1,1.73)));
    c = max(c, p.x);
    return c;
}


vec4 HexCoords(vec2 uv) {
	vec2 r = vec2(1, 1.73);
    vec2 h = r*.5;
    vec2 a = mod(uv, r)-h;
    vec2 b = mod(uv-h, r)-h;
    vec2 gv = dot(a, a) < dot(b,b) ? a : b;
    float x = atan(gv.x, gv.y);
    float y = .5-HexDist(gv);
    vec2 id = uv-gv;
    return vec4(x, y, id.x,id.y);
}
#define pi acos(-1.)
#define tau (2.*pi)
vec3 pCoords = vec3(0);

float tCoords = 0.;
vec3 cCoords = vec3(1.);

#define coolPal(a,b) pal(0.5,0.6,vec3(0.97 + 0.2*sin(p.z*0.04),3.4 + sin(b)*0.2,0.8),0.4 + (a),3.3 + (b))
#define tunnW 0.8


float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); }

vec2 map(vec3 p){
	vec2 d = vec2(10e6);
    
    //p.xy *= rot(0. + sin(p.z)*0.24 + cos(p.y)*0.13);
    
	p -= path(p.z);
    p.xy *= rot(0. + sin(p.z)*0.24 + cos(p.y)*0.13);
    //p.xy *= rot(0.6);
	//float wave = pow(abs(sin(p.z*0.2 + 0.5*((smooth_basstime*0.457)+(TIME*0.125)))), 40.);
    //float wave = pow(abs(sin(p.z*0.2 + smoothTimeC*0.125 + 0.17)), 40.);
    float wave = pow(abs(sin(p.z*0.2 + smoothTimeC*0.125 + 0.17)), 20.)*(0.5+syn_HighLevel);
    
    pCoords = vec3(atan(p.y,p.x)/tau, length(p.xy), p.z);
    
    float tunn = -length(p.xy + normalize(p.xy)*wave*0.15) + tunnW;
    
    //tunn = max(tunn, - max(abs(p.x)- 0.86,abs(p.y) - 0.3 ));
    
    d = dmin(d, vec2(tunn, 2.));
    
    p.xy *= rot(0.5 + sin(p.z*0.2 + p.y*0.4)*1.);
    
    vec3 pCoordsB = pCoords;
    pCoordsB.y -= 0.4;
    pCoordsB.x = pmod(pCoordsB.x, 0.05)/0.05;
    pCoordsB.z = pmod(pCoordsB.z, 1.)/1.;
    
    pCoordsB = abs(pCoordsB) - 0.5;
    float dCubeGlow = max(pCoordsB.x, max(pCoordsB.y,pCoordsB.z));
    
    
    glow += max(exp(-dCubeGlow*10.)*0.024*pow(wave,0.14)*pal(0.5,0.5,vec3(1.6,2.9,0.9),0.5,0.5)*(0.7+.5*syn_Intensity), 0.);
    
    float attA = pow(abs(sin(sin(p.z + (smoothTimeB*0.2)) + p.z*0.7 + sin(p.y*1.5)) ), 10.);
    float attB = pow(abs(sin(smoothTimeB*0.2 + sin(p.z + (smoothTimeB*0.2)) + p.z*0.4 + sin(p.y*1.5))), 10.);
    glow += attB*smoothstep(0.,1.,length(p.xy)*0.8)*max(exp(-dCubeGlow*12.)*0.02*pal(0.8,0.6,vec3(0.3,0.6,0.5),0.5,0.4)*2.1, 0.)*(0.5+syn_HighLevel+syn_MidHighLevel*0.5);

    
    vec3 c = max(coolPal(01.3, 1.5), 0.);
    
    vec3 q = p;
    float pipe;

    
    
    q = p;
    q.xy -= 0.4;
    pipe = length(q.xy) - 0.076;
    d = dmin(d, vec2(pipe,10.));
    q.xy += 0.9;
    pipe = length(q.xy) - 0.076;
    d = dmin(d, vec2(pipe,10.));
    
    
    
    float dtt = 10e7;
    float dsqP = 10e6;
    q = p;
    q.z -= 0.125*smoothTime;

    q.xy -= vec2(-0.6,0.2);
    vec3 z = abs(vec3(q.x,q.y,pmod(q.z, 0.4))) - vec3(0.1,0.04,0.2);
    z.xy *= rot(0.125*pi);
    float rc = max(z.y,max(z.x, z.z));
    q = abs(q) - 0.07;
    dtt = rc + 0.02 ;
    pipe = max(q.x,q.y);
    pipe = max(pipe, -rc);
    dsqP = min(dsqP, pipe);
    
    q = p;
    q.z -= 0.125*smoothTime;

    q.xy += vec2(-0.4,0.4);
    z = abs(vec3(q.x,q.y,pmod(q.z, 0.4))) - vec3(0.1,0.04,0.2);
    z.xy *= rot(0.125*pi);
    rc = max(z.y,max(z.x, z.z));
    q = abs(q) - 0.07;
    dtt = min(dtt,rc + 0.02 );
    pipe = max(q.x,q.y);
    pipe = max(pipe, -rc);
    dsqP = min(dsqP, pipe);
    
    glow += exp(-dtt*(10. - pow(wave,0.4)*10. + sin(p.z)*1.))*max(coolPal(0.8, 1.8), 0.)*1.*pow(abs(sin(p.z + (.5*smoothTimeB))), 10.);
    glow += exp(-dtt*(50. - pow(wave,0.4)*10. + sin(p.z)*1.))*max(coolPal(0.8, 1.8), 0.)*(0.9+syn_HighLevel);
    
    d = dmin(d, vec2(dtt,6.));
    
    
    q = p;
    q.z = pmod(q.z, 0.2);
    
    d = dmin(d, vec2(dsqP,6.));
    
    float net;
    
    pCoordsB = pCoords;
    pCoordsB.y -= 0.7;
    pCoordsB.x = pmod(pCoordsB.x,1./8.)/(1./8.);
    pCoordsB.z = pmod(pCoordsB.z,1.)/(1.);
    float ww = 0.01;
    net = length(pCoordsB.xy) - ww;
    net = min(net,length(pCoordsB.yz) - ww);
    
    glow += exp(-net*(20. - pow(wave,0.4)*10. + sin(p.z)*4.))*max(coolPal(0.7, 1.8), 0.);
    d.x *= 0.7;
	return d;
}

vec2 march(vec3 ro, vec3 rd,inout vec3 p,inout float t,inout bool hit){
	vec2 d = vec2(10e6);
	t = 0.; hit = false; p = ro;
    for(int i = 0; i < 96; i++){
    	d = map(p);
        glow += exp(-max(d.x, 0.)*20.)*0.01;
        if(d.x < 0.001){
        	hit = true;
            break;
        }
    	t += d.x;
        p = ro + rd*t;
    }
    
    
    return d;
}
vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
	vec3 dir = normalize(lookAt - ro);
	vec3 right = cross(vec3(0,1,0), dir);
	vec3 up = cross(dir, right);
    float fov = 1.;
	return normalize(dir + right*uv.x*fov + up*uv.y*fov);
}
vec3 getNormal(vec3 p){
	vec2 t = vec2(0.001, 0.);
    return normalize(map(p).x - vec3(
    	map(p-t.xyy).x,
    	map(p-t.yxy).x,
    	map(p-t.yyx).x
    ));
}


#define mx (10.*iMouse.x/RENDERSIZE.x)
#define my (10.*iMouse.y/RENDERSIZE.x)
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    vec3 col = vec3(0);
	//uv *= 1. + dot(uv,uv)*1.5;
    
    vec3 ro = vec3(0);
    ro.z += (bass_time);
    ro.z += mx;
        
    ro += path(ro.z);
    
    float wave = pow(abs(sin(ro.z*0.2 + smoothTimeC*0.125 + 0.17)), 40.);
    
    vec3 lookAt = vec3(0,0,ro.z);

    lookAt.xy += _uvc*PI*Fov*PI;
    
    lookAt.z += 2.;
    lookAt += path(lookAt.z); 
	vec3 rd = getRd(ro, lookAt, uv*(1. + wave*0.016));
        rd.xy = _rotate(rd.xy, Rotate*PI);

    ro -= rd*texture(image30,(uv)*256.).x*0.2; // remove banding from glow
    
    
    //rd.xy *= rot(sin(script_time*0.12 + sin(script_time*0.14)*0.125 )*0.1);
    rd.yz = _rotate(rd.yz, lookXY.y*PI+PI*_uvc.y*Flip*Fov);
    rd.xz = _rotate(rd.xz, lookXY.x*PI);

    float t; vec3 p; bool hit;
    
    vec2 d = march(ro, rd, p, t, hit);
     
    
    if(hit){

    	vec3 n = getNormal(p);       
        if(d.y == 2.){
            vec4 hc = HexCoords(vec2(8.*pCoords.x, pCoords.z*1.)*1.);
            vec4 hcc = hc;
            float dHex = 10e6;
            
            float md = 0.1 ;
            
            float formula = hc.y + 0.8*sin(p.z*0.03 + (smoothTimeC*0.0125))+0.4*cos(p.z+0.0125*smoothTimeC);
            formula += hc.w*0.1;
            float id = floor(formula/md);
            hc.y = (pmod(formula, md))/md;
            
            
            dHex = min(dHex,abs(hc.y) - 0.1);
            
            vec3 c = coolPal(0. + 0.8, sin(id*0.7)*0.19);
            c *= max(0.,sin(id*0.9 + sin(hc.w)*9.));
            
			c -= hcc.y*0.7;
            c = max(c, 0.);
            
            col += c;
            
            float atten = sin(p.x);
            
        }
        if(d.y == 8.){
            float dd;
            
            cCoords = pmod(abs(cCoords),0.2)/0.2;
            
            cCoords = cCoords - 0.1;
        }
        if(d.y == 6.){
            float dd;
            
            col += pow(1. - max(dot(n, -rd), 0.), 5.)*0.3;
            
            
        }
        if(d.y >= 10.){
        	vec3 c = pal(0.6,0.1,vec3(01.7,2.4,3.3), 0.2,3.3 + sin(p.z*0.5)*0.5);
            
            c = max(c, 0.);
            float a = pmod(p.z + (smoothTimeC*0.0125), 1.)/1.;
            col += pow(1. - max(dot(n, -rd), 0.), 5.)*0.1;

            col += smoothstep(0.03,0.,abs(a) - 0.14)*c;
        }
        
    }
    
    col += pow(glow,vec3(1./0.45))*0.001;
    
    col = max(col, 0.);
    //col = mix(col, vec3(0.4,0.5,0.9)*1.5,pow(smoothstep(0.,1.,t*0.08), 3.9));
 
	col = mix(col, vec3(0.7,0.3,0.2)*3.5,pow(smoothstep(0.,1.,t*0.08), 6.9));   
    col *= 0.896;
    
    
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


// Fork of "Day 76" by jeyko. https://shadertoy.com/view/WtyXzt
// 2020-03-05 21:05:41

// hex function from BigWIngs - check out his tutorial on youtube

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    
    
    //float m = pow(abs(sin(p.z*0.03)),10.);

    // Radial blur
    float steps = 20.;
    float scale = 0.00 + pow(dot(uvn,uvn),1.)*0.04;
    float chromAb = dot(uvn*1.5,uvn*1.5)*10.;
    vec2 offs = vec2(0) + texture(BuffA, uv + smoothTimeC*.125).xz*0.0;
    vec4 radial = vec4(0);
    for(float i = 0.; i < steps; i++){
    
        scale *= 0.99;
        vec2 target = uv + offs;
        offs -= normalize(uvn)*scale/steps;
    	radial.x += texture(BuffA, target + chromAb*1./RENDERSIZE.xy).x;
    	radial.y += texture(BuffA, target).y;
    	radial.z += texture(BuffA, target - chromAb*1./RENDERSIZE.xy).z;
    }
    
    
    
    radial /= steps;
    
    fragColor = texture(BuffA, uv)*0.6 + radial*0.4; 
    
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.);
    fragColor *= 0.7;
    
    fragColor = pow(fragColor, vec4(0.9));
    
    fragColor.b *= 0.5;
    fragColor = smoothstep(0.,1.,fragColor);
    //fragColor = pow(fragColor, vec4(0.45));
	
    float duvuv = dot(uvn,uvn);
    
    fragColor = pow(fragColor, vec4(0.4545 + 0.87*clamp(duvuv*0.4, 0., 0.50)));
    
    fragColor *= 1. - duvuv*1.4;
    
    //fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.5);
    
    
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