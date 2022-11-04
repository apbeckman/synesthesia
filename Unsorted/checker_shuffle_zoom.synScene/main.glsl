

// From FabriceNeyret2 (262 chars)
//*
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 u = _xy;

    vec2 m, R = RENDERSIZE.xy,a;
    O*=0.;
    for (float x,e, i = 0.; i < 9.; m.x*m.y < 0. ? O += .11 : O )
        x = .004 * i++ - smoothTime*0.2,
        e = fract(x+x), 
        a = abs( m = mod( 9.*(u-R/2.)/R.y * ++e, 2.) - 1. ) -.5,
        fract(x) > .5 ? m = m.yx : m,
        m *= sign(m.y),
        a.y * a.x < 0. ? m.x -= sign(a.y)/e : e,
        m = fract(m) - .5;
        return O;
}
/**/

// From Xor (266 chars)
/*
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 u = _xy;

    vec2 m, R = RENDERSIZE.xy,a;
    O=vec4(0);
    for (float x,e, i = 0.; i < 9.; m.x*m.y < 0. ? O += .11 : O )
        x = .004 * i++ - TIME,
        e = fract(x+x), 
        a = abs( m = mod( 9.*(u-R/2.)/R.y * ++e, 2.) - 1. ) -.5,
        fract(x) > .5 ? m = m.yx : m,
        m *= sign(m.y),
        a.y * a.x < 0. ? m.x -= sign(a.y)/e : e,
        m = fract(m) - .5;
}
/**/

// From FabriceNeyret2 (267 chars)
/*
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 u = _xy;

    vec2 m, R = RENDERSIZE.xy,a;
    O=vec4(0);
    for (float x,e, i = 0.; i < 9.;  m.x*m.y < 0. ? O += .11 : O )
        x = .004 * i++ - TIME,
        e = fract(x+x), 
        a = abs( m = mod( 9.*(u-R/2.)/R.y * ++e, 2.) - 1. ) -.5,
        fract(x) > .5 ? m = m.yx : m,
        m.y > 0. ? m = -m : m,
        a.y * a.x < 0. ? m.x += sign(a.y)/e : e,
        m = fract(m) - .5;
}
/**/

// From Xor (291 chars)
/*
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 u = _xy;

    vec2 m, R = RENDERSIZE.xy;
    u -= R * .5;
    O = vec4(0);

    for (float x, y, e, a, A = 9., s = A / R.y, i = 0.; i < A; O += vec4(m.x > .5 ^^ m.y > .5) / A)
    
        x = fract(.004 * i++ - TIME) * s * 2.,
        y = mod(x, s),
        e = s + y,
        a = y / e,
        m = mod(u * e, 2.),
        x > s ? m = m.yx : m,
        m.y > 1. ? m = 2. - m : m,
        m.y > .5 ^^ abs(m.x - 1.) < .5 ? m.x -= m.y > .5 ? -a : a : a,
        m = fract(m);
}
/**/

// My initial golfing attempt (315 chars)

/*
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 u = _xy;

    vec3 R = RENDERSIZE;
    u -= R.xy * .5;
    O = vec4(0);

    for (float A = 9., s = A / R.y, i = 0.; i < A; ++i)
    {
        float x = fract(.004 * i - TIME) * (s + s),
            y = mod(x, s),
            e = s * (1. + y / s),
            a = y / e;
            
        vec2 m = mod(u * e, 2.);
        if (x > s) m = m.yx;
        if (m.y > 1.) m = 2. - m;
        if (m.y > .5) a = -a;
        if (m.y > .5 ^^ abs(m.x - 1.) < .5) m.x -= a;
        m = fract(m);
        
        O += vec4(m.x > .5 ^^ m.y > .5) / A;
    }
}
/**/


// Original: 492 chars

/*
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 u = _xy;

    u -= RENDERSIZE.xy * 0.5;
    float size = 10.0 / RENDERSIZE.y;
    
    float AA = 16.;
    O = vec4(0);
    for (float i = 0.; i < AA; ++i)
    {
        float factor = -(TIME - 0.016667 * (i / AA)) * 2.0;
        float mm = mod(factor * size, size * 2.);
        float mm0 = mod(mm, size);
        float scale = size * (1.0 + mm0 / size);
        float animate = mm0 / scale;
        vec2 m = mod(u * scale, 2.0);
        if (mm > size) m = m.yx;
        if (m.y > 1.0) m = 2.0 - m;
        if (m.y < 0.5 && m.x > 0.5 && m.x < 1.5) m.x += animate;
        if (m.y > 0.5 && (m.x > 1.5 || m.x < 0.5)) m.x -= animate;
        m = mod(m, 1.0);
        O += vec4(m.x > 0.5 ^^ m.y > 0.5) / float(AA);
    }
	return O; 
 } 

/**/


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}