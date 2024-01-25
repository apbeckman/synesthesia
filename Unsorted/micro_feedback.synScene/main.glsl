			//******** BuffA Code Begins ********
#define T(x) texture(BuffA, fract((x)/RENDERSIZE.xy))

vec4 renderPassA() {   
    vec4 c = vec4(0.0);
    vec2 u = _xy;
    //c=1./u.yyyx;

    c=u.yyyx/1e4;///iTime;
    //for(float t=1.4; t<1e2; t+=t)
    //    c += (c.gbar-c)/3.+T(u-c.wz*t);
    for(float t=.6; t<4e2; t+=t)
    	c += c.gbar/4.-c*.3+T(u-c.wz*t);
            c += _textureMedia(_uv);
    u -= _uvc;
	c = mix(T(u-0.5*c.xy), cos(c), .07);
	return c;
}

vec4 renderMainImage() {
    vec4 c = vec4(0.0);
    vec2 u = _xy;
    c = .5+.5*texelFetch(BuffA, ivec2(u),0);
	return c;
}


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}