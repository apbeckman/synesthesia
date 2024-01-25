			//******** Common Code Begins ********
float dist(vec2 p0, vec2 pf) {
    return sqrt((pf.x - p0.x) * (pf.x - p0.x) + (pf.y - p0.y) * (pf.y - p0.y));
}
float smin(float a, float b, float k) {
    float h = max(k - abs(a-b), 0.) / k;
    return min(a, b) - h*h*h*k*1./6.;
}

float d2 = dist(RENDERSIZE.xy * 0.5, _xy.xy) * 0.004;
float shock = (shock_back - shock_forward*(1.0+syn_Intensity));
vec2 rotate(vec2 coords, float angle) {
    float sin_factor = sin(angle);
    float cos_factor = cos(angle);
    coords = vec2((coords.x - 0.5), coords.y - 0.5) * mat2(cos_factor, sin_factor, -sin_factor, cos_factor);
    coords += 0.5;
    return coords;
}
float twist = spin_right - spin_left;
vec2 uvc_noise = vec2(_noise(_uv.x), _noise(_uv.y));
float zm_amt = 1.0 - (1.0+low*0.2)*_scale(Zoom, -.0025, .0025);
#define R RENDERSIZE.xy
vec2 mirror(vec2 u) {
    if(u.x > 1.)
        u.x = 1. - fract(u.x);
    if(u.x < 0.)
        u.x = fract(u.x);
    if(u.y > 1.)
        u.x = 1. - fract(u.y);
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

vec4 renderPassA() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;
    Q = C(U);
    vec4 n = D(U + vec2(0, 1)), e = D(U + vec2(1, 0)), s = D(U - vec2(0, 1)), w = D(U - vec2(1, 0)), m = 0.25 * (n + e + s + w);
    float d = 0.25 * (n.y - s.y + e.x - w.x);
    float c = 0.25 * (n.x - s.x - e.y + w.y);
    float zm_amt = 1.0 + _scale(Zoom, -0.01,0.04);
    // Q.z = m.z * (.999*zm_amt) - mix(d, c, length(U - 0.5 * R) / R.y);
    Q.z = m.z * .999 - mix(d, c, length(U - 0.5 * R) / R.y);


    Q.w = d;
    if(FRAMECOUNT <= 1)
        Q = vec4(sin(U.x) * cos(U.y));
    return Q;
} 

			//******** BuffB Code Begins ********

vec4 renderPassB() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;
    U += _uvc * Stretch * PI * (1.40 + 0.5 * low);
    U -= _uvc * PI * (Zoom - shock * 1.5) * (1.40 + low);
    U += _uvc * PI * (Fisheye + shock * PI) * (1.0 + 0.5 * low) * d2;
    U = rotate(U, 0.00 + (twist * 0.005) * d2);
    U +=uvc_noise*warp;
    U += vec2(dot(cos(_uvc.x), sin(_uvc.y))) * Pull * PI * 2;
    U += mix(_rotate(vec2(_noise(abs(mix(_uvc, _uv, 0.5) + 0.1 * TIME) * PI)), _noise(smoothTimeC * 0.025)) * (warp + stir * PI * 2), vec2(_noise(abs(_uv + 0.2 * TIME) * PI)) * (warp + stir * PI), noise_rot);

    Q = A(U);
    vec4 n = A(U + vec2(0, 1)), e = A(U + vec2(1, 0)), s = A(U - vec2(0, 1)), w = A(U - vec2(1, 0)), m = 0.25 * (n + e + s + w);
    Q.xy = 0.25 * vec2(e.z - w.z, n.z - s.z);
    
    if(length(Q.xy) > 0.)
        Q.xy = mix(Q.xy, normalize(Q.xy), .2);
    return Q;
} 

			//******** BuffC Code Begins ********

vec4 renderPassC() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    Q = A(U);
    vec4 n = B(U + vec2(0, 1)), e = B(U + vec2(1, 0)), s = B(U - vec2(0, 1)), w = B(U - vec2(1, 0)), m = 0.25 * (n + e + s + w);

    float d = 0.25 * (n.y - s.y + e.x - w.x);
    float c = 0.25 * (n.x - s.x - e.y + w.y);

    Q.z = m.z * .999 - mix(d, c, .2);

    if(FRAMECOUNT <= 1)
        Q = vec4(sin(U.x) * cos(U.y));
    return Q;
} 

			//******** BuffD Code Begins ********

vec4 renderPassD() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;

    vec2 c = 0.5 * R;
    // if(_mouse.z > 0.)
    //     c = _mouse.xy;
    c.xy *= 1.0 + ((moveXY.xy )*(1.0+low));

    U -= c;
    U *= .99*zm_amt;
    U += c;
    Q = C(U);
    vec4 n = C(U + vec2(0, 1)), e = C(U + vec2(1, 0)), s = C(U - vec2(0, 1)), w = C(U - vec2(1, 0)), m = 0.25 * (n + e + s + w);
    Q.xy = 0.25 * vec2(e.z - w.z, n.z - s.z);
    if(length(Q.xy) > 0.)
        Q.xy = mix(Q.xy, normalize(Q.xy), .2);

    return Q;
}

mat2 r(float a) {
    vec2 e = vec2(cos(a), sin(a));
    return mat2(e.x, -e.y, e.y, e.x);
}
vec4 renderMainImage() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;


    vec3 p = vec3(0.5 * R, .62 * R.y), d = normalize(vec3((U - 0.5 * R) / R.y, -1)), o = vec3(0.5, .1, .5) * R.xyy;
    // if(_mouse.z > 0.)
    //     o.xy = _mouse.xy;
    mat2 m = r(.44);
    p.y -= .19 * R.y;
    float m_x = smin(d.x*-1, d.x, 0.125);
    float m_y = smin(d.y*-1, d.y, 0.125);
    vec2 rd_mirrored = vec2(m_x, m_y);
    d.xy = mix(d.xy, rd_mirrored, vec2(mirror_x, mirror_y));

    d.yz *= m;
    p.yz *= m;
    for(int i = 0; i < 35; i++) {
        p += .2 * d * (p.z - 4. * A(p.xy).z);
    }
    d = normalize(o - p);
    float z = A(p.xy).z;
    vec3 n = normalize(vec3(B(p.xy).xy, -.25));
    vec3 q = d;
    p += .1 * d;
    for(int i = 0; i < 35; i++) {
        p += .5 * d * min(p.z - 4. * A(p.xy).z, length(p - o) - 1.);
    }
    Q = (exp(-length(p - o) + 1.)) * (cos(-.1 * TIME + .1 * z + .45 * vec4(1, 2, 3, 4))) * .5 * (dot(reflect(n, d), q) - dot(n, d));
    Q *= Q;
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