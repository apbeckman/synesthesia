

			//******** Common Code Begins ********

vec2 R;
float T;
#define Main void mainImage(out vec4 Q, in vec2 U)
#define ei(a) mat2(cos(a),-sin(a),sin(a),cos(a))
float ln (vec4 p, vec4 a, vec4 b) {
    float l = clamp(dot(p-a,b-a)/dot(b-a,b-a),0.,1.);
    return mix(.8,1.,l)*length(p-a-(b-a)*l);
}
vec4 map (vec4 u) {
    vec4 Q =vec4(0);
    vec2 U = _xy;

    u.xz *= ei(T);
    u.xy *= ei(1.5);
    float d = 1e9;
    vec4 c;
    float sg = 1e9;
    float l = .2;
    u.y = abs(u.y);
    u.y+=.2;
    for (float i = 0.; i < 12.; i++)
    {
        sg = ln(u,vec4(0),vec4(0,l,0,0))/l;
        d = min(d,sg);
        u.y -= l;
        u.xz *= ei(2.);
        u.xz = abs(u.xz);
        u.xy *= ei(.3);
        l *= .8;
            
    }
    return vec4(1.)*smoothstep(.1,0.,abs(d-.5));;
}


// Fork of "Microscopy 101" by wyatt. https://shadertoy.com/view/NlKSzG
// 2021-12-31 17:42:39


vec4 renderMainImage() 
{
    R = RENDERSIZE.xy;
    T = TIME*0.1;
    vec4 Q =vec4(0);
    vec2 U = _xy;
    for (float i = -100.; i < 100.;i++){
        vec4 u = vec4(U,i/150.,0);
        u.xy = 1.*(u.xy-.5*R.xy)/R.y;
    
        Q += map(u)/30.;
    }
	return Q; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}