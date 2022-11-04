#define 	TWPI  	6.283185307179586  	// two pi, 2*pi

vec4 texMirror(sampler2D samplerIn, vec2 uvIn){
    if (mod(uvIn.x, 2.0) > 1.0){
        uvIn.x = 1.0-uvIn.x;
    }
    if (mod(uvIn.y, 2.0) > 1.0){
        uvIn.y = 1.0-uvIn.y;
    }
    return texture(samplerIn, mod(uvIn, 1.0));
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    vec2 v = (_uv.xy + RENDERSIZE.y) - 0.5;

    float th = v.y * PI,
        ph = v.x * TWPI;
    vec3 sp = vec3( sin(ph) * cos(th), sin(th), cos(ph) * cos(th) );
    vec3 pos = vec3( PI, 0, 0);

    pos *= ((mix(0.25+(1.0-zoom)*4., 0.25+(1.0-zoom)*2., perspective)) * 0.1);
    //pos += (zoom * 0.18 * perspective);
    sp = mix(sp, normalize(vec3(_uvc, 1.0)), 1.0-perspective);

    sp.yz = _rotate(sp.yz, look_xy.y*PI);
    sp.xy = _rotate(sp.xy, look_xy.x*PI);

    vec2 nUv = vec2(dot(pos, sp.zxy), dot(pos.yzx, sp.zxy));
    nUv -= 0.5;
    nUv.x += velocity.x;
    nUv.y += velocity.y;
    nUv += (position);
    nUv.y += mix(0.0, 1.0, mirror);

    if (_exists(syn_UserImage)) {
        fragColor = _contrast(_invertImage(mirror > 0.0 ? texMirror(syn_UserImage, nUv) : texture(syn_UserImage, mod(nUv, 1.0))),_Media_Contrast);
    } else {
        fragColor = texture(image1, mod(nUv * 6.0, 1.0));
    }

    return fragColor;
}
vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}
