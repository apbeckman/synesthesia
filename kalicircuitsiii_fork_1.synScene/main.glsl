

mat2 rot(float a) {
	float s=sin(a), c=cos(a);
    return mat2(c,s,-s,c);
}


vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}


vec3 fractal(vec2 p) {
    float o=100.,l=o;
    p*=rot(-.3);
    p.x*=1.+p.y*.7;
    p*=.5+sin(TIME*.1)*.3;
    p+=TIME*.02;
    p=fract(p);
    for (int i=0; i<10; i++) {
        p*=rot(radians(90.));
        p.y=abs(p.y-.25);
        p=p/clamp(abs(p.x*p.y),0.,3.)-1.;
        o=min(o,abs(p.y-1.)-.2)+fract(p.x+TIME*.3)*.5;
        l=min(l,length(max(vec2(0.),abs(p)-.5)));
        
    }
	o=exp(-5.*o);
	l=smoothstep(.1,.11,l);
    return hsv2rgb(vec3(o*.5+.6,.8,o+.1))+l*vec3(.4,.3,.2);
}


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy;
    uv-=.5;
    uv.x*=RENDERSIZE.x/RENDERSIZE.y;
    vec3 col=vec3(0.);
    float aa=4.;
    vec2 eps=1./RENDERSIZE.xy/aa;
    for (float i=-aa; i<aa; i++){
        for (float j=-aa; j<aa; j++){
            col+=fractal(uv+vec2(i,j)*eps);
        }
    }
	col/=pow(aa*2.,2.);
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}