vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


// -----------------------------------------------------
// wheelofconf. by nabr
// https://www.shadertoy.com/view/3dXfRH
// License Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
// https://creativecommons.org/licenses/by-nc/4.0/
// -----------------------------------------------------


float pch(vec2);
#define t TIME
float f(vec3 p)
{
    float f = 0., 
          z = .05 * t + p.z, 
          v = 1. - smoothstep(0., 2., length(p));
    for(float i = 10.; i > 0.; --i){
        float c = p.y * sin(i / 1.57 + z), 
              s = p.x * cos(i / 1.57 + z), 
              b = v * mod(1.57 * v, .636);
        f += min(v, .0125 / abs(b +(fract(.05*t)>.5?2.:1.) * (c + s)));
    }
    return (p + normalize(smoothstep(0., f, vec3(.15, 1, .5)))).z;
}
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 U = _xy;

    vec2 R = RENDERSIZE.xy, e = vec2(.001, 0);
    U = (U - .5 * R) / R.y;
    vec3 c, p, N, o = vec3(0, 0, .5+(sin(.2*t)*.2+.2)), d = vec3(U, -1),
                  l = -normalize(normalize(vec3(.7 * sin(.5 * t), cos(.5 * t), 1)) + d),
                  dn = vec3(0, 0, 6.f), D = dot(dn, o) / dot(dn, d) * d - o;
    float h, lr = exp2(.05 * dot(l, l));
    for(int i = 10; i > 0; --i)h += f(p = (D * h));
    N = normalize(f(p) - vec3(f(p - e.xyy), f(p - e.yxy), f(p - e.yyx)));
    c = mix(vec3(.2, .6, .95), max(dot(N, l), .05) * lr * vec3(.8), fract(.05 * t + .5 + h) + .5);
    c += vec3(.1, .2, .2) * max(0., .8 - dot(reflect(l, N), -d)) * fract(.5 * h);
    O = vec4(vec3(.4)* pch(U)+clamp(c, 0., 1.), 1.);
    return O;
}
// char by FabriceNeyret2
#define char(_p,_C)(((_p.x<0.||_p.x>1.||_p.y<0.||_p.y>1.))?vec4(0,0,0,1e5):textureGrad(iChannel3,_p/16.+fract(vec2(_C,15-_C/16)/16.),dFdx(_p/16.),dFdy(_p/16.))).x
int msg[]=int[](84,104,101,32,119,111,114,108,100,32,119,105,108,108,32,115,116,105,108,108,32,98,101,32,116,117,114,110,105,110,103,32,119,104,101,110,32,121,111,117,39,114,101,32,103,111,110,101,89,101,97,104,44,32,119,104,101,110,32,121,111,117,39,114,101,32,103,111,110,101,45,32,79,122,122,121);
float pch(vec2 u){
    vec2 tp=24.*vec2(u.x+.25,u.y-.3);
    float p; 
    for(int i=(msg.length()-1);i>=0;--i){
        p+= char(tp,msg[i]);
        switch(i){
        case 70:tp.x-=2.;tp.y-=.75;break;
        case 48:tp.x-=24.;tp.y-=.75;
        }
        tp.x+=.5;
    };return p*max(0., 1.-t*.1); 
	return O; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}