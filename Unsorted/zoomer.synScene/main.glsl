

			//******** Common Code Begins ********
float dist(vec2 p0, vec2 pf) {
    return sqrt((pf.x - p0.x) * (pf.x - p0.x) + (pf.y - p0.y) * (pf.y - p0.y));
}

float d = dist(RENDERSIZE.xy * 0.5, _xy.xy) * (_mouse.x / RENDERSIZE.x + 0.1) * 0.005;
float d2 = dist(RENDERSIZE.xy * 0.5, _xy.xy) * 0.00125;

#define R RENDERSIZE.xy
#define A(U) texture(BuffA,(U)/R)
#define B(U) texture(BuffB,(U)/R)
#define C(U) texture(BuffC,(U)/R)
#define D(U) texture(BuffD,(U)/R)
#define Main void mainImage(out vec4 Q,vec2 U)
#define box for(int x=-1;x<=1;x++)for(int y=-1;y<=1;y++)
#define _a 1e-4
#define ei(a) mat2(cos(a),-sin(a),sin(a),cos(a))

			//******** BuffA Code Begins ********

vec4 renderPassA() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;
    U -= _uvc * Zoom * (1.0 + low * syn_Intensity)*PI;
    U -= _uvc * Stretch;
    U += Drift;
    U -= _uvc*d2*fisheye*PI;

    Q = A(U);
    vec4 Qb = B(U);
    vec2 du = vec2(0), m = vec2(0);
    box if(abs(x) != abs(y)) {
        vec2 u = vec2(x, y);
        vec4 a = A(U + u);
        vec4 b = B(U + u);
        du -= .25 * mat2(b) * u;
        m += .25 * a.xy;
    }
    Q.xy = mix(Q.xy, m, 1.);
    Q.xy = mix(Q.xy, du, .33);
    if(length(Q.xy) > 0.)
        Q.xy = mix(Q.xy, normalize(Q.xy), .3);
    if(FRAMECOUNT <= 1)
        Q = vec4(sin(.04 * U.y), cos(.04 * U.x), 0, 0);
    return Q;
} 

			//******** BuffB Code Begins ********

vec4 renderPassB() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;
    U -= _uvc * Zoom;
    U += Drift;
    U -= _uvc * Stretch;

    Q = B(U);
    vec4 Qa = A(U);
    vec4 du = vec4(0), m = vec4(0);
    box if(abs(x) != abs(y)) {
        vec2 u = vec2(x, y);
        vec4 a = A(U + u);
        vec4 b = B(U + u);
        du += vec4(u.x * a.xy, u.y * a.xy);
        m += .25 * b;
    }
    Q = mix(Q, m, 1.);
    Q = mix(Q, du, .1);
    if(FRAMECOUNT <= 1)
        Q = vec4(sin(.01 * U.y), -cos(.01 * U.x), 0, 1);
    return Q;
} 

			//******** BuffC Code Begins ********

vec4 renderPassC() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;
    U -= _uvc * Zoom;
    U += Drift;
    U -= _uvc * Stretch;

    U -= _uvc * Zoom;
    U -= 0.5 * R;
    U *= (1. - _a);
    U += 0.5 * R;
    Q = A(U)*(1.0+syn_Level);
    return Q;
} 

			//******** BuffD Code Begins ********

vec4 renderPassD() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;
    U -= _uvc * Zoom;
    U += Drift;
    U -= _uvc * Stretch;

    U -= 0.5 * R;
    U *= (1. - _a);
    U += 0.5 * R;
    Q = A(U);
    return Q;
}

vec4 renderMainImage() {
    vec4 Q = vec4(0.0);
    vec2 U = _xy;
    U -= _uvc * Zoom;
    U += Drift;
    U -= _uvc * Stretch;
    U -= _uvc*d2*fisheye*PI;
    vec4 a = A(U + _uvc), b = B(U + _uvc);
    float t = b.x * b.y + b.w * b.z;
    Q = a.xzyw + b + t * vec4(10 * (1.0 + pow(syn_HighLevel * 0.5 + syn_Intensity * 0.35 + syn_MidHighLevel * 0.35, 2.)));
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