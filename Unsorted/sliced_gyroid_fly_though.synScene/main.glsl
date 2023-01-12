

/** 
    License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
    
    Sliced Gyroid Fly Though
    11/16/22 | @byt3_m3chanic

    music for mood only - Greenhouse Gasses | by Cevin Key

    Original post here about sliced SDF's. Kind of updated
    version (extra loops to help backside of sdf's)
    
    @smjtyazdi https://twitter.com/smjtyazdi/status/1484828390104485896

*/

#define R   RENDERSIZE
#define M   _mouse

#define PI  3.1415926

mat2 rot (float a) { return mat2(cos(a),sin(a),-sin(a),cos(a)); }
float hash21(vec2 a){ return fract(sin(dot(a,vec2(22.34,35.34)))*483434.);}

//globals
mat2 trot;
float glow=0.,tglow=0.;

const float ofx = 1.25;

float sdGry(vec3 p, float s, float t, float b) {
    p *= s;
    return abs(dot(sin(p*ofx),cos(p.zxy))-b)/(s*ofx)-t;
}

vec2 map(vec3 p) {
    vec2 res =vec2(1e5,1.);
    
    vec3 q = p;
    float tt = .5+.5*sin(TIME*.2);
    
    float k = 8./dot(p,p);
    float mul = 1./k;
    p*=k;
    
    p.xy*= trot;
    p.yz += vec2(-1.25, bass_time*2.);
    
    float d = .365+Sliced, mf = 1e5;
    float mm = 1.5+.5*sin(high_time*0.2);
    
    for(float j=-2.;j<2.;j++){
        vec3 nf =p;
        nf.xy += _uvc;
        nf.z=round(nf.z/d+j)*d;
        float ids = mod(nf.z,2.);
        
        float fd= sdGry(nf, .725, .05, .65);
        
        nf.z=clamp(p.z,nf.z-d/4.,nf.z+d/4.);
        fd=length(vec2(max(.0,fd), nf.z-p.z));

        float idx = mod(nf.z,5.);
        
        if(idx<mm+.025&&idx>mm) tglow+=.002/(.015+fd*fd);
        if(ids<mm+.075&&ids>mm) glow+=.002/(.015+fd*fd);
   
        mf=min(mf,fd);
    }
    if(mf<res.x) res=vec2(mf,3.);

    res.x*= mul/1.35;
    return res;
}

vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 F = _xy;


    trot=rot(spin_time*.072+Spin*PI);
    
    vec2 uv = (2.* F.xy-R.xy)/max(R.x,R.y);
    uv += _uvc*PI*FOV;
    vec3 ro = vec3(0,0,.5);
    vec3 rd = normalize(vec3(uv, -1));
    rd.xz = _rotate(rd.xz, LookXY.x*PI);
    rd.yz = _rotate(rd.yz, LookXY.y*PI);
    float x = M.xy == vec2(0) || M.z < 0. ? 0. : -(M.y/R.y * .5 - .25) * PI;
    float y = M.xy == vec2(0) || M.z < 0. ? 0. : -(M.x/R.x * 2. - 1.) * PI;

    mat2 rx = rot(x),ry = rot(y);
    
    ro.yz *= rx;ro.xz *= ry;
    rd.yz *= rx;rd.xz *= ry;
    
    vec3 C = vec3(0);
    vec3 p = ro;

    float d=0.,m=0.;
    for(int i=0;i<164;i++){
        vec2 t = map(p);
        d += i<64? t.x*.35:t.x*.75;
        m  = t.y;  
        p = ro + rd * d;
        if(t.x<d*1e-4||d>45.) break;
    } 
    
    float sp = .2+.2*sin(uv.x*4.1*PI+bass_time);
    //vec3 fog = mix(vec3(0.6043,0.253,0.255)*0.125,vec3(0.235,0.02,0.9000)*0.25,clamp((uv.y+.5-sp),0.,1.))*pow(1.0+syn_HighLevel*(1.0+syn_Intensity), 2.0);
    vec3 fogCol = vec3(0.143,0.1253,0.255);
    vec3 fog = mix(fogCol*0.125,fogCol*0.25,clamp((uv.y+.5-sp),0.,1.))*pow(1.0+syn_HighLevel*(1.0+syn_Intensity), 2.0);
    fog.xy += _uvc*PI*0.01;

    C = mix(C,vec3(glow,glow*.65,glow*.15),clamp(glow,.0,.6));
    C = mix(C,vec3(tglow*.15,tglow*.465,tglow),clamp(tglow,.0,.8));
    
    C = mix(C,fog, 1.-exp((Fog*-20.)*d*d*d));
    
    C = clamp(C,vec3(0),vec3(1));
    O = vec4(pow(C, vec3(.4545)),1);
	return O; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}