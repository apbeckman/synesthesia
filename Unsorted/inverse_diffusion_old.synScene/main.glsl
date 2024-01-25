//******** Common Code Begins ********
float dist(vec2 p0, vec2 pf) {
    return sqrt((pf.x - p0.x) * (pf.x - p0.x) + (pf.y - p0.y) * (pf.y - p0.y));
}
float d2 = dist(RENDERSIZE.xy * 0.5, _xy.xy) * 0.004;
float shock = shock_back - shock_forward;
vec2 rotate(vec2 coords, float angle) {
    float sin_factor = sin(angle);
    float cos_factor = cos(angle);
    coords = vec2((coords.x - 0.5), coords.y - 0.5) * mat2(cos_factor, sin_factor, -sin_factor, cos_factor);
    coords += 0.5;
    return coords;
}
#define R RENDERSIZE.xy
vec2 mirror(vec2 u) {
    if(u.x > 1.)
        u.x = 1. - fract(u.x);
    if(u.x < 0.)
        u.x = fract(u.x);
    if(u.y > 1.)
        u.x = 11. - fract(u.y);
    if(u.y < 0.)
        u.x = fract(u.y);
    return u;
}
#define A(U) texture(BuffA,mirror((U)/R))
#define B(U) texture(BuffB,mirror((U)/R))
#define C(U) texture(BuffC,mirror((U)/R))
#define Ca(U) texture(syn_UserImage,mirror((U)/R))*(0.5+low)
#define D(U) texture(BuffD,mirror((U)/R))

//******** BuffA Code Begins ********
// vec2 Mirror = vec2(x_flip, y_flip);

vec4 renderPassA() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    U -= _uvc * PI * Zoom * (1.0 + low);

    Q = C(U);
    vec4 n = D(U + vec2(0, 1 + Height)), e = D(U + vec2(1 + Height + distortion, 0)), s = D(U - vec2(0, 1 + Height + distortion)), w = D(U - vec2(1, 0)), m = 0.25 * (n + e + s + w);

    float d = 0.25 * (n.y - s.y + e.x - w.x);
    float c = 0.25 * (n.x - s.x - e.y + w.y);

    Q.z = m.z * .9999 - mix(d, c, length(U - 0.5 * R) / R.y);
    Q.w = d;

    if(FRAMECOUNT <= 1 || Reset != 0.)
        Q = vec4(sin(U.x) * cos(U.y));
    return Q;
}

//******** BuffB Code Begins ********
float twist = spin_right - spin_left;

vec4 renderPassB() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;
    vec2 uvc_noise = vec2(_noise(_uvc.x), _noise(_uvc.y));
    // U += _rotate(vec2(_noise(abs(_uv + 0.1 * TIME) * PI)), _noise(smoothTimeC*0.05)) * warp;
    U += mix(_rotate(vec2(_noise(abs(mix(_uvc, _uv, 0.5) + 0.1 * TIME) * PI)), _noise(smoothTimeC * 0.095)) * (warp + stir * PI * 2), vec2(_noise(abs(_uv + 0.3 * TIME) * PI)) * (warp + stir * PI), noise_rot);
    U += vec2(dot(cos(_uvc.x), sin(_uvc.y))) * Pull * PI * 2;
    U += _uvc * Stretch * PI * (1.40 + 0.5 * low);
    U -= _uvc * PI * (Zoom - shock * 1.5) * (1.40 + 0.5 * low);
    U += _uvc * PI * (Fisheye + shock * PI) * (1.0 + 0.5 * low) * d2;
    U = rotate(U, 0.00 + (twist * 0.005) * d2);

    Q = A(U);
    vec4 n = A(U + vec2(0, 1)), e = A(U + vec2(1, 0)), s = A(U - vec2(0, 1)), w = A(U - vec2(1, 0)), m = 0.25 * (n + e + s + w);
    Q.xy = 0.25 * vec2(e.z - w.z, n.z - s.z);

    if(length(Q.xy) > 0.)
        Q.xy = mix((Q.xy), normalize(Q.xy), .25);
    return Q;
}

//******** BuffC Code Begins ********

vec4 renderPassC() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    Q = A(U);
    vec4 n = B(U + vec2(0, 1 + distortion)), e = B(U + vec2(1 + distortion, 0)), s = B(U - vec2(0, 1 + distortion)), w = B(U - vec2(1 + distortion, 0)), m = 0.25 * (n + e + s + w);

    float d = 0.25 * (n.y - s.y + e.x - w.x);
    float c = 0.25 * (n.x - s.x - e.y + w.y);

    Q.z = m.z * .9999 - mix(d, c, .01);
    if(FRAMECOUNT <= 1 || Reset != 0.)
        Q = vec4(sin(U.x) * cos(U.y));
    // Q += -0.075*edgeFeedback+mix(vec4(0), abs(texture(syn_FinalPass, _uv)*0.075), edgeFeedback);

    return Q;
}

//******** BuffD Code Begins ********

vec4 renderPassD() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;
    //U += ((moveXY.xy*(1.0+0.5*low)));

    vec2 c = (R);
    c.xy *= 1.0 + ((moveXY.xy * (1.0 + 0.5 * low)));

    // c.xy = mix(abs(c.xy), c.xy, 1.0 - Mirror);

    c *= 1.0 - vec2(_noise(TIME * 0.7), _noise(TIME * 0.7));
    // c += _noise(c);
    // if (_mouse.z>0.) c = _mouse.xy;
    // U.xy = mix(abs(U.xy), U.xy, 1.0 - Mirror);

    U -= c;
    //    U += _uvc*PI;

    //    U *= (.99325*(1.-growthFactor*Zoom*0.00675)-Zoom*0.01);
    U *= (.999);

    U += c;
    Q = C(U);

    vec4 n = C(U + vec2(0, 1 + Height)), e = C(U + vec2(1 + Height, 0)), s = C(U - vec2(0, 1 + Height)), w = C(U - vec2(1 + Height, 0)), m = 0.25 * (n + e + s + w);
    Q.xy = 0.25 * vec2(e.z - w.z, n.z - s.z);
    //  Q.r += _pixelate(texture(syn_FinalPass, _uv).r, .125)*12*edgeFeedback;

    if(length(Q.xy) > 0.)
        Q.xy = mix(Q.xy, normalize(Q.xy), .1);

    return Q;
}
vec4 feedback() {
    vec4 Q = _edgeDetectSobel(syn_FinalPass, _uv);
    vec2 U = _xy;
    return Q;
}
float circle = length(_uvc);
mat2 r(float a) {
    vec2 e = vec2(cos(a), sin(a));
    return mat2(e.x, -e.y, e.y, e.x);
}
vec4 renderMainImage() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;
    vec3 p = vec3(0.5 * R + _uvc, .62 * R.y + _uvc.y * PI), d = normalize(vec3((U + _uvc - 0.5 * R) / R.y, -1)), o = vec3(0.5, .1, .5) * R.xyy;
    o.xy += vec2(cos(TIME), sin(TIME));
    // if (_mouse.z > 0.)
        // o.xy = _mouse.xy;
    mat2 m = r(.44);
    p.xy += vec2(dot(cos(_uv.x + _uvc.x), sin(_uv.y + _uvc.y))) * Pull * PI;
    p.xy += vec2(_noise((_uv)));
    p.y -= .19 * R.y;
    d.yz *= m;
    p.yz *= m;
    for(int i = 0; i < 35; i++) {
        p += (.125) * d * (p.z - (5.) * A(p.xy).z);
    }

    float highs = (0.35 * pow(syn_HighLevel * 0.35 + syn_MidHighLevel * 0.35 + syn_Intensity * .2, 2));
    d = normalize(o - p);
    float z = A(p.xy).z;
    vec3 n = normalize(vec3(B(p.xy).xy, -1.4)); //.4

    vec3 q = d;

    p += .1 * d;
    for(int i = 0; i < 24; i++) {
        p += .5 * d * min(p.z - 5. * A(p.xy).z, length(p - o) - 1.);
        Q += smoothstep(0.1, 0.4, circle) * vec3(_pulse(circle, fract(TIME + 3.33 * i), 0.1));
    }
    vec4 baseCol = vec4(1, 2, 3, 4);
    //Q = (exp(-length(p-o)+1.)*(1.)*(cos(-.015*smoothTimeC*.2+(.15*(1.)) * z + .5*vec4(1,2,3,4)*-1.))*.5*(dot(reflect(n,d),q)-dot(n,d)))*(1.0+pow(syn_HighLevel, 2)*0.1);
    // vec4 baseCol_adjusted = _hueRotate(baseCol, hue_rot);
    Q = (exp(-length(p - o) + 1.) * (cos(-.05 * smoothTimeB + (.175 * (1.5)) * z + .5 * baseCol * -1.)) * .5 * (dot(reflect(n, d), q) - dot(n, d)));
    // Q.rgb = _palette(0.925, vec3(Q.rgb), vec3(0.5, 0.5, 0.5	), vec3(1.0, 1.0, 1.0), vec3(0.00, 0.33, 0.67));
    Q = _contrast((_grayscale(Q, greyscale * 1.)), .975);
    Q = _brightness(Q, 1.0 + (syn_MidHighLevel * 0.05 + syn_HighLevel * 0.04) * flash);
    // Q = _hueRotate(Q, hue_rot*180);
    Q *= Q;
    // Q.rgb += texture(syn_FinalPass, _uv).rgb*0.1256*edgeFeedback;
    return Q;
}

vec4 renderMain() {
    if(PASSINDEX == 0) {
        return renderPassA();
    }
    if(PASSINDEX == 1) {
        return renderPassB();
    }
    if(PASSINDEX == 2) {
        return renderPassC();
    }
    if(PASSINDEX == 3) {
        return renderPassD();
    }
    if(PASSINDEX == 4) {
        return renderMainImage();
    }
}
