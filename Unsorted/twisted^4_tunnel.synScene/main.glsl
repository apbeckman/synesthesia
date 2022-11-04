vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


#define MDIST 60.0
#define STEPS 178.0
#define pi 3.1415926535
#define rot(a) mat2(cos(a),sin(a),-sin(a),cos(a))
#define pmod(p,x) (mod(p,x)-0.5*(x))


vec3 hsv(vec3 c){
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}
//My poor mans version of Javad Taba's helix function
vec3 spiral(vec3 p, float R){
   
    p.xz*=rot(p.y/R);
    vec2 s = sign(p.xz);
    p.xz=abs(p.xz)-R*0.5;
    
    p.xz*=rot(smoothTimeC*0.5*pi/3.);
    
    float poy = p.y;
    p.y=0.;
    //This is supposed to correct for the distortion that happens
    //when you twist a shape by rotating it over an axis.
    //In my head it should work a lot better than it does, but
    //it definitely helps a little bit so I guess it's better than nothing
    p.yz*=rot(mix(0.,pi/4.,1./(R*0.5+1.5)))*-sign(s.x*s.y);
    p.y=poy;
    return p;
}


vec2 map(vec3 p){
    float t = smoothTime*0.25;
    //p.y+=sin(-p.z*0.1)*2.;
    p.y-=p.z*p.z*0.008;
    
    p.zy*=rot(pi/2.);
    vec3 po = p;

    p.y-=t*pi*4.339;
    
    vec2 a = vec2(1);
    vec2 b = vec2(2);
    p.xz*=rot(-0.05*(smoothTimeC/6.));
    p = spiral(p,6.6);
    p = spiral(p,2.);
    p = spiral(p,1.);
    p = spiral(p,0.4);
    //there are some small artifacts but you dont notice them ;) 
    
    
    //p = spiral(p,0.1);
    //vec2 d = abs(p.xz);
    //a.x = max(d.x,d.y)-1.0;
    a.x = length(p.xz)-0.1;
    //a.x = max((abs(po.y)-7.),a.x);
    a.x*=0.6;
    return vec2(a);
}
vec3 norm(vec3 p){
    vec2 e = vec2(0.005,0);
    return normalize(map(p).x-vec3(
    map(p-e.xyy).x,
    map(p-e.yxy).x,
    map(p-e.yyx).x));
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    vec3 col = vec3(0);
    vec3 ro = vec3(0,0.,-1.)*1.5;
    if(iMouse.z>0.){
    ro.yz*=rot(3.0*(iMouse.y/RENDERSIZE.y-0.5));
    ro.zx*=rot(-7.0*(iMouse.x/RENDERSIZE.x-0.5));
    }
    vec3 lk = vec3(0,0,0);
    vec3 f = normalize(lk-ro);
    vec3 r = normalize(cross(vec3(0,1,0),f));
    vec3 rd = normalize(f*(0.9)+uv.x*r+uv.y*cross(f,r));  
    vec3 p = ro;
    float dO = 0.;
    bool hit = false;
    vec2 d= vec2(0);
    for(float i = 0.; i<STEPS; i++){
        p = ro+rd*dO;
        d = map(p);
        dO+=d.x;

        if(abs(d.x)<0.005||i>STEPS-1.5){
            hit = true;
            break;
        }
        if(dO>MDIST){
            dO = MDIST;
            break;
        }
    }
    if(hit)
    {
        vec3 ld = normalize(-vec3(p.x,p.y,p.z-5.));
      
        //sss from nusan
        float sss=0.01;
        for(float i=1.; i<20.; ++i){
            float dist = i*0.35;
            sss += smoothstep(0.,1.,map(p+ld*dist).x/dist)*0.18*1.25;
        }
        for(float i=1.; i<5.; ++i){
            float dist = i*0.7;
            sss += smoothstep(0.,1.,map(p-ld*dist).x/dist)*0.25;
        }
        vec3 al = vec3(0.204,0.267,0.373);
        vec3 n = norm(p);
        vec3 r = reflect(rd,n);
        float diff = max(0.,dot(n,ld));
        float amb = dot(n,ld)*0.45+0.55;
        float spec = pow(max(0.,dot(r,ld)),40.0);
        #define AO(a,n,p) smoothstep(-a,a,map(p+n*a).x)
        float ao = AO(.3,n,p)*AO(.5,n,p)*AO(.9,n,p);

        col = al*
        mix(vec3(0.169,0.000,0.169),vec3(0.984,0.996,0.804),mix(amb,diff,0.75))
        +spec*0.3;
        col+=sss*hsv(vec3(0.76,0.9,1.35));
        col*=mix(ao,1.,0.5);
        col = pow(col,vec3(0.7));
    }
    vec3 bg = mix(vec3(0.094,0.000,0.200),vec3(0.600,0.000,0.600),length(rd.xy)-0.65);
    
    col = mix(col,bg,pow(dO/MDIST,2.5));
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}