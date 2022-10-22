

			//******** Common Code Begins ********


//#define R RENDERSIZE.xy

#define T(u) texture(iChannel0,(u)/R)
#define T1(u) texture(iChannel1,(u)/R)
#define T2(u) texture(iChannel2,(u)/R)
#define T3(u) texture(iChannel3,(u)/R)


#define Neighbors vec4 n = T(U + vec2(0,1)); vec4 s = T(U - vec2(0,1)); vec4 e = T(U + vec2(1,0)); vec4 w = T(U - vec2(1,0)); vec4 m = 0.25*(n + w + e + s);




//MK2022 = Anaglyph Sphere Waves =
//This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 

#define MAX_STEPS 200
#define MAX_DIS 400.
#define MIN_DIS 10.

#define SURF_DIS .01
#define SURF_MUL 100.
#define SURF_EXP 2.

#define EYE_DIS 1.0
#define FIL_COR 1.05

#define TIME TIME*.6

mat2 Rot(float a)  
{
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

float Dist(vec3 p) 
{
    //cam rotation
    float q = min(smoothTimeC/200., 5.);
    float c = 0.2 / (1.+q/3.);
    p.yx *= Rot(.7854*smoothstep(0.,c+c,-cos((smoothTime*0.1)/(23.+q))));
    p.zx *= Rot(.7854*smoothstep(-c,c+c,-cos((smoothTime*0.1)/(11.+q))));
    p.zy *= Rot(.7854*smoothstep(c,c+c,-cos((smoothTime*0.1)/(45.+q))));

    //distance frequency
    float f = 1./(400.+sin(smoothTimeC/23.)*200.+sin(smoothTimeC/29.)*100.);
    float dis = length(p)*f;
    
    //mirrorbox
    vec3 size = vec3(10.);
    p = mod(abs(-mod(p,4.*size)+2.*size),4.*size);
    
    //sphere
    float d = length(p - sin(smoothTimeC*.5+dis*6.28)*size.x-size) - 
              10.*(sin(smoothTimeC*1.+dis*6.28)*.5+.8)*pow(dis/f/MAX_DIS, 1.2) *
              (1.+smoothstep(0.9,1.,sin(smoothTimeC/2.+dis/f/200.))); 

    return d;
}

float RTM(vec3 ro, vec3 rd) 
{
    float sum = 0.;
	float s = 1.;
    float d = MIN_DIS;
    const float a = 1. / float(MAX_STEPS); 
    
    for(int i = 0; i < MAX_STEPS; i++) 
    {    
        float sd = (SURF_DIS * (pow(d/MAX_DIS, SURF_EXP)*SURF_MUL+1.));
        if (s < sd || d > MAX_DIS) break;
        
        s = Dist(ro + rd*d);
        s = max(abs(s), 2.*sd);
        d += s * .8;
        sum += a;
    }
    
    return smoothstep(0., 1., sum * (1.-exp(-d*d)));
}

vec3 R(vec2 uv, vec3 p, vec3 l, float z)
{
    vec3 f = normalize(l-p),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f,r),
        c = p+f*z,
        i = c + uv.x*r + uv.y*u,
        d = normalize(i-p);
    return d;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-.5*RENDERSIZE.xy) / RENDERSIZE.y;
	
    //interaction
    vec2 m = _mouse.xy / RENDERSIZE.xy;
    if (length(m) == 0.) m = vec2(.5);
    m.x = clamp(m.x, .1, .9);
    m.y = clamp(m.y, .25, .75);
    
    //set cam
    vec3 ro = vec3( 0., 0., -15.);
    ro.yz *= Rot(-m.y*3.14159 + 1.5708);
    ro.xz *= Rot(-m.x*6.28318);
    vec3 rd = R(uv, ro, vec3(0., 0., 0.), .7 + 0.4*smoothstep(.6,.9,sin(smoothTime/16.)));
    ro += normalize(cross(rd, vec3(0., -1., 0.))) * EYE_DIS * -.5;
    
    //left eye
    float colL = RTM(ro, rd);
    
    //eye distance
    ro += normalize(cross(rd, vec3(0., -1., 0.))) * EYE_DIS;
    
    //right eye
    float colR = RTM(ro, rd);
    
    //filter correction
    colL *= FIL_COR; 
    colR /= FIL_COR;
    
    //contrast
    vec3 colS = smoothstep(vec3(0.), vec3(1.), vec3(colL, colR, colR));
    
    //Vignette
    colS *= smoothstep(2.,-2./5., dot(uv,uv)); 
    
    fragColor = vec4(colS, 1.);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}