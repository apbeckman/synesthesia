mat3 rot (vec3 s) {
    float 	sa = sin(s.x),
        ca = cos(s.x),
        sb = sin(s.y),
        cb = cos(s.y),
        sc = sin(s.z),
        cc = cos(s.z);
    return mat3 (
        vec3(cb*cc, -cb*sc, sb),
        vec3(sa*sb*cc+ca*sc, -sa*sb*sc+ca*cc, -sa*cb),
        vec3(-ca*sb*cc+sa*sc, ca*sb*sc+sa*cc, ca*cb)
    );
}
mat3 mm;
vec4 light;
vec3 grp;
float myAVR;
float ui;
float map (vec3 p) {
    float d = length(p-light.xyz)-light.w;
    d = min(d,max(30.-p.z,0.));
    grp=p*0.01;
    
    p = mm*p;
    for (int i = 0; i < 5; i++) {
        float t = abs(p.y);
        p.y = p.x;
        p.x = t; //zoom??
        p = mm*(p-.1*multiply);
    }
    float q = 5.;
    
    float a = 1.5*((atan(p.x,p.z)))+0.05*(ui);
    mat2 mn = mat2(sin(a),cos(a),-cos(a),sin(a));
    vec2 u = mn*vec2(length(p.xz)-2.-expand,p.y);
    float d1 = d;
    vec2 w = max(abs(u)-vec2(0.1,0.5),0.);
    d = min(d,length(w)-0.1-0.01*(sin(50.*(dot(u,u))+sin(20.*atan(p.x,p.z)))));
    if (d1 != d) grp = vec3(u,atan(p.x,p.z));
    return d;
}
vec3 norm (vec3 p) {
    vec2 e = vec2 (.001,0.);
    return normalize(vec3(
        map(p+e.xyy) - map(p-e.xyy),
        map(p+e.yxy) - map(p-e.yxy),
        map(p+e.yyx) - map(p-e.yyx)
    ));
}
vec3 dive (vec3 p, vec3 d) {
    for (int i = 0; i < 35; i++) {
        p += d*map(p);
    }
    return p;
}
vec3 color (vec3 no, vec3 p) {
    return vec3(.8+0.2*(sin((30.)*p*p.zxy+no)*0.5+0.5))*0.9;
}
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
    	vec2 fragCoord = _xy;
    vec2 v = (-RENDERSIZE.xy + 2.0*fragCoord.xy)/RENDERSIZE.y;
    ui = (rate*script_time*0.5)*20.;
    vec3 r = vec3(1-(1.5*moveXY.x),0-moveXY.y,-15.);
    light = vec4(10.*sin(0.005*ui),2,-23,1);
    vec3 d = normalize(vec3(v,3.+zoom));
    mm = rot((0.005*ui)+vec3(1,2,3));
    vec3 p = dive(r,d);
    d = normalize(light.xyz-p);
    vec3 no = norm(p);
    vec3 col = color(no,grp);
    vec3 bounce = dive(p+0.01*d,d);
    col *= dot(no, normalize(light.xyz-(p)));
    if (length(bounce-light.xyz) > light.w+0.1) col *= 0.2;
    if (length (p-r)>4e2) col *= 0.;
    fragColor = vec4(col,1);
    return fragColor;
}

vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}