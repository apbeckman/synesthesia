

			//******** Common Code Begins ********



mat2 rot(float angle){
    return mat2(cos(angle), -sin(angle), sin(angle), cos(angle));
}
float hash13(vec3 p3){
    p3 = fract((p3)*0.1031);
    p3 += dot(p3, p3.yzx  + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}



float valueNoiseo(vec3 uv,float pw){
    vec3 id = floor(uv);
    vec3 fd = fract(uv);
    fd = smoothstep(0.,1., fd);
    
    fd = pow(fd,vec3(pw));
    
    float ibl = hash13(id + vec3(0,-1,0));
    float ibr = hash13(id + vec3(1,-1,0));
    float itl = hash13(id + vec3(0));
    float itr = hash13(id + vec3(1,0,0));
    
    
    float jbl = hash13(id + vec3(0,-1,1));
    float jbr = hash13(id + vec3(1,-1,1));
    float jtl = hash13(id + vec3(0,0, 1));
    float jtr = hash13(id + vec3(1,0, 1));
    
    
    float ibot = mix(ibl, ibr, fd.x); 
    float iup = mix(itl, itr, fd.x);
    float jbot = mix(jbl, jbr, fd.x);
    float jup = mix(jtl, jtr, fd.x);
    
    float i = mix(ibot, iup, fd.y);
    float j = mix(jbot, jup, fd.y);
    
    return mix(i, j, fd.z); 
}



			//******** BuffA Code Begins ********

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    fragColor = vec4(0.0,0.0,1.0,1.0);
	return fragColor; 
 } 



// CHEAP NOISE FROM NIMITZ https://www.shadertoy.com/view/XtS3DD
// RENDERING TECHNIQUE LEARNED FROM IQ'S "CLOUDS"
#define TIME (TIME + _mouse.x/RENDERSIZE.x)

float valueNoise(in vec3 p,float pw)
{
    vec3 ip = floor(p);
    vec3 fp = fract(p);
	fp = fp*fp*(3.0-2.0*fp);
	vec2 tap = (ip.xy+vec2(37.0,17.0)*ip.z) + fp.xy;
	vec2 rz = textureLod( image30, (tap+0.5)/256.0, 0.0 ).yx;
	return mix( rz.x, rz.y, fp.z );
}


float fbm(vec3 p){

    vec3 op = p;
    float n = 0.;
    p *= 0.5;
    p.y += TIME*0.5;
    float fa = valueNoise(p,1.); 
    
    p.y += fa*(1. + sin(op.z + TIME*0.2)*0.5);
    p.x += TIME*.325;
    p.y -= TIME*0.25;
    
    float fb = valueNoise(p*2.,1.);
    
    
    p.x += fb*.2 + TIME*0.1;
    p.z += fb*1.1 + TIME*0.02;
    float fd = valueNoise(p*9.8,1.);
    
    float fc = valueNoise(p*4.2,1.);
    
    

    n = fa*1. + fb*0.55 + fc*0.244 + fd*0.146*(0.6 - fa) /*fa*/;// + valueNoise(op*5.8,1.)*0.0;
    
    
    //n*= 1.1;
    
    //n = smoothstep(0.,1.,n);
    //n = pow(n,1. + fa*2.);
    
    
    //n *= pow(valueNoise(op + 24.,4.),0.4)*3.;
    
    
    n*= 0.6;
    
    n = smoothstep(0.,1.,n);
    n = pow(n,1. + fa*1.5);
    
    n *= pow(valueNoise(op*0.3 + 1.,1.),0.4)*1.;
    n *= 2.4;
    return clamp( (1. - op.y*0.1 + abs(op.x)*0.4 - 3.0 + 2.*n) , 0.0, 1.0 );

	return n;
}

//vec3 ro;

float map(vec3 p){
    float d = 10e5;
    vec3 op = p;
    
    p.xy *= rot(p.z*0.14);
    
    
    p.z += TIME;
    p.y += sin(TIME)*0.2 + 0.;
    //p.xy *= rot(p.z*0.1);
    
    p.xy *= rot(sin(p.z*0.2 + TIME*0.2)*0.4);
    
    d = fbm(p);
    
    
    float dopop = dot(op,op);
    
    d *= smoothstep(0.,1., (length(p.xy)*0.6 - 2.4)*smoothstep(1.,0.,
        dopop*.03) + smoothstep(0.,1., dopop*.01)
        );
    
    
    return d;
}
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    vec3 col = vec3(0.12,0.14,0.14)*1.6;
    col = vec3(0.14,0.14,0.12)*1.6;


    col += sin(uv.xyx*0.2 + TIME*0.4)*0.1;
    
    
    col += sin(uv.xyx*0.2 )*0.1;
    
    vec3 lightDir = normalize(vec3(-1.,1.,2.));

    mat2 lightRot = rot(TIME*0.5); 
    
    lightDir.xy *= lightRot;
    vec3 ro = vec3(0);
    vec3 rd = normalize(vec3(uv,0.9));
    
    rd.yz *= rot(sin(TIME*0.25)*0.2);
    
    vec3 p = ro;
    
    
    float stSz = 0.1;
    float dTotal = 0.;
    vec3 accum = vec3(0);
    
    ro += 0.2*rd*hash13(vec3(uv*220.,1.+ TIME*4.));
    
    vec3 lightCol;
    
    
    
    float steps = 200.;
    float i = 0.;
    float t = 0.;
    for(; i < steps; i++){
        float dCurr = map(p);
        float dOther = map(p + lightDir*.15 );
        
        lightCol = vec3(0.4 + sin(TIME*0.25 + p.z*0.0)*0.2,0.2,0.1);
        vec3 colDiff = vec3(0.115,0.4,0.4)*0. + col*0.7 + 32.7*lightCol*clamp((dCurr - dOther*1.)*0.6,0.,1.); 
        vec3 colAbsorption = mix(vec3(0.47 + sin(TIME + p.z)*0. - 0.2,0.5,.65)*0.5, vec3(0.4,0.2,0.9)*0.2, smoothstep(0.,1.,dTotal)); 
        vec3 colFringe = clamp( (1. - dCurr*2.)*lightCol*vec3(0.8,0.9,.9 )*0.2, 0., 1.);
        
        if(dTotal > 1. || t > 20222.){
            break;
        }
        
        
        float fade = pow(smoothstep(0.,1.,(p.z - ro.z)*0.03),0.75);
        
        vec3 colCurr = mix( (colDiff)*colAbsorption + colFringe*1.9 , col, fade);
        dCurr *=  1. -fade;
        accum += dCurr*colCurr*(1. - dTotal)*stSz;
        dTotal += dCurr*stSz*(1. - dTotal);
        t += max(stSz,i/steps*2.*stSz);
        //p += rd*t;
        p = ro + rd*t;
    }


    //lightCol = vec3(0.4 + sin(TIME*0.45)*0.2,0.2,0.1);
        
    
    col = mix(col,accum*1.,pow(dTotal, 1.) );
    col += smoothstep(1.,0.,length( (uv + vec2(0.7,-0.4)*lightRot) ) + 0.1)*lightCol*(0.2 + 0.*accum*accum*dTotal);
    
    
    col = mix(col,smoothstep(0.,1.,col*1.7),0.7);
    col = mix(col,smoothstep(0.,1.,col*1.6),1.);
    
    col *= smoothstep(1.,0.3,dot(uv,uv)*0.6 + 0.2);
    col = pow(max(col,0.),vec3(0.45454));
    
    fragColor = vec4(col,1.0);
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