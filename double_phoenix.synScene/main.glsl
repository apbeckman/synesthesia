

			//******** Common Code Begins ********


#define R RENDERSIZE.xy

#define T(u) texture(iChannel0,(u)/R)
#define T1(u) texture(iChannel1,(u)/R)
#define T2(u) texture(iChannel2,(u)/R)
#define T3(u) texture(iChannel3,(u)/R)


#define Neighbors vec4 n = T(U + vec2(0,1)); vec4 s = T(U - vec2(0,1)); vec4 e = T(U + vec2(1,0)); vec4 w = T(U - vec2(1,0)); vec4 m = 0.25*(n + w + e + s);




//A 'Double Phoenix' fractal
//Julia is Zn+1 = Zn^2 + c
//Phoenix is Zn+1 = Zn^2 + c + a * Zn-1
//I added a term of b * Zn-2, just for fun
//a, b, and c real-only for symmetry
//It looks nice
//Using Morgan's heatmap gradient
//https://www.shadertoy.com/view/4dsSzr

vec3 heatmapGradient(float t) {
    t = 2.*abs(.5-fract(t));
	return clamp((pow(t, 1.5) * 0.8 + 0.2) * vec3(smoothstep(0.0, 0.35, t) + t * 0.5, smoothstep(0.5, 1.0, t), max(1.0 - t * 1.7, t * 7.0 - 6.0)), 0.0, 1.0);
}

float phoenix(vec2 p, vec3 c){
	vec2 p1 = vec2(0), p2=p1;
    for(int i=0;i<200;i++){
     	vec2 tmp = p, sq = p*p;
        if(sq.x + sq.y > 65536.*2.0)
            return float(i) - log2(log2(sq.x+sq.y)*.07625);
        p = vec2(sq.x-sq.y+c.x,(2.+test)*p.x*p.y) + c.y*p1 + c.z*p2;
        p2 = p1;
        p1 = tmp;
    }
    return -1.;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = 3.*(fragCoord.xy -.5 * RENDERSIZE.xy) / RENDERSIZE.y;
    float t = phoenix(uv.yx, sin(smoothTime*0.2*vec3(.13,.17,.19)));
    vec3 col = t > 0. ? heatmapGradient(log(t) + .1*smoothTime) : vec3(.05);
	fragColor = vec4(col ,1.0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}