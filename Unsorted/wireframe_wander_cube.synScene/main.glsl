

			//******** BuffA Code Begins ********

// Fork of "cube wireframe thingy" by Tater. https://shadertoy.com/view/sddGR8
// 2021-08-17 03:45:12

// Fork of "Lost Wander Cube" by Tater. https://shadertoy.com/view/Nst3RH
// 2021-08-17 02:45:45

//Add more steps to destroy your PC and make it look cooler
#define STEPS 60.0

#define MDIST 100.0
#define pi PI
#define rot(a) mat2(cos(a),sin(a),-sin(a),cos(a))
#define pmod(p, x) (mod(p,x)-0.5*(x))

float box(vec3 p, vec3 s){
    vec3 d = abs(p)-s;
    return max(d.x,max(d.y,d.z));
}
float frame(vec3 p, vec3 s, float e){
    vec2 h = vec2(e,0);
    float a = box(p,s);
    float b = box(p,s*1.01-h.xxy);
    a = max(-b,a);
    b = box(p,s*1.01-h.xyx);
    a = max(-b,a);
    b = box(p,s*1.01-h.yxx);
    a = max(-b,a);
    return a;
}
float timeRemap (float t,float s1, float s2, float c){
    return 0.5*(s2-s1)*(t-asin(cos(t*pi)/sqrt(c*c+1.0))/pi)+s1*t;  
}
void mo(inout vec2 p){
  if(p.y>p.x) p = p.yx;
}
vec3 kifs;
vec2 map(vec3 p){
    vec3 po2 = p;
    
    p.xz*=rot(TIME*0.8);
    p.xy*=rot(TIME*0.4);
    vec3 po = p;
    float t = (TIME+50.0)*0.7;
    t = timeRemap(t*1.3, 0.0, 2.3, 0.1);

    for(float i = 0.0; i< 9.0; i++){
        p = abs(p)-2.0*i*kifs;
        p.xz*=rot(pi/2.0);
        mo(p.xy);
        mo(p.zy);
       // p.x-=sign(p.y)*sin(t);
       
    }
    
    //Inner Cubes
    p = pmod(p,2.2+sin(smoothTimeC*0.5)*0.5+0.5*cos(smoothTimeC*0.5));
    vec2 a = vec2(box(p,vec3(0.5)),1.0);
    a.x = abs(a.x)-0.2;
    a.x = abs(a.x)-0.1;
    vec2 b = vec2(box(p,vec3(0.45)),2.0);
    
    a = (a.x<b.x)?a:b;
    
    p = po;
    p.xy*=rot(pi/4.0);
    
    //Boundry Cut Cube
    vec3 cube = vec3(4,4,4)*vec3(1.2+0.5*sin(t),1.2+0.5*cos(t),1.2+0.5*sin(t));
    a.x = max(box(p,cube),a.x);
    b.x = max(box(p,cube),b.x);
   
    //Outer Frame
    b= vec2(frame(p,cube+0.15,0.45),3.0);
    a = (a.x<b.x)?a:b;
    return a;
}
vec3 norm(vec3 p, float s) {
  vec2 off=vec2(s,0);
  return normalize(
   vec3(map(p+off.xyy).x, map(p+off.yxy).x, map(p+off.yyx).x)
  -vec3(map(p-off.xyy).x, map(p-off.yxy).x, map(p-off.yyx).x));
}
void render( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (fragCoord-0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    vec3 col = vec3(0);
    vec3 ro = vec3(0,2,-19);
    vec3 lk = vec3(0,0,0);
    vec3 f = normalize(lk-ro);
    vec3 r = normalize(cross(vec3(0,1,0),f));
    vec3 rd = f*1.0+uv.x*r+uv.y*cross(f,r);
    float t = (TIME+50.0)*0.7;
    t = timeRemap(t*1.3, 0.0, 2.3, 0.1);
    kifs = abs(vec3(asin(sin(t*0.15)),0.2*asin(sin(t*0.22)),0.3*asin(sin(t*0.38))));
     float s=1.0;
    
    float rayDist,shad;
    vec3 p = ro;
    vec2 d;
    bool hit = false;
    vec3 al;
    
    //Raymarch loop
    for(float i =0.0; i<STEPS; i++){
        p = ro+rd*rayDist;
        d = map(p);
        if(abs(d.x)<0.01){
            float edge = 0.04*clamp(800.0/RENDERSIZE.x,1.0,2.0);
            float edgeAmount = length(norm(p, 0.015)-norm(p, edge));
            
            //Adjusting line transparency based on steps
            col += smoothstep(0.0,0.1, edgeAmount)*0.045*(128.0/min(STEPS,200.0));
            s*=-1.0;
            d.x = 0.01*s;
            if(!hit){
            if(d.y == 1.0)al = vec3(1.000,0.302,0.404);
            if(d.y == 2.0)al = vec3(0.016,1.000,1.000);
            if(d.y == 3.0)al = vec3(1.000,0.647,0.325);
                shad = i/STEPS;
                hit = true;
            }
        }
        if(rayDist>MDIST) break;
        rayDist+=d.x*s*0.85;
       
    }
    vec3 wire = col;
    wire = clamp(wire,0.0,1.0)*al;
    shad = pow(1.0-shad,2.0);
    if(hit)col = vec3(shad)*sqrt(al);
    col = wire;
    fragColor = vec4(col,1.0);
}
#define AA 1.0
vec4 renderPassA() {
	vec4 O = vec4(0.0);
	vec2 C = _xy;

    float px=1./AA,i,j;vec4 cl2,cl;
    if(AA==1.){render(cl,C);O=cl;return O;}
    for(i=0.;i<AA +min(TIME,0.0);i++){for(j=0.;j<AA;j++){
    vec2 C2 = vec2(C.x+px*i,C.y+px*j);
    render(cl2,C2);cl+=cl2;
    }}cl/=AA*AA;O=cl;
	return O; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    fragColor = texture(BuffA,fragCoord/RENDERSIZE.xy,0.);

    for(float i = 0.; i < 7.; i+=0.24){
        fragColor += texture(BuffA, fragCoord/RENDERSIZE.xy, i)/18.;
    }
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