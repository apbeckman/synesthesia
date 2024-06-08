vec3 iMouse = vec3(xyPaint1 * RENDERSIZE, paint1_on);
float dist(vec2 p0, vec2 pf){return sqrt((pf.x-p0.x)*(pf.x-p0.x)+(pf.y-p0.y)*(pf.y-p0.y));}
float d = dist(RENDERSIZE.xy*0.5,_xy.xy)*(_mouse.x/RENDERSIZE.x+0.1)*0.005;
float d2 = dist(RENDERSIZE.xy*0.5,_xy.xy)*0.00125;

#define _G0 0.25
#define _G1 0.125
#define _G2 0.0625
#define W0 -3.0
#define W1 0.5
#define TIMESTEP 0.1
float ADVECT_DIST = mix(2.0, 2.0 - syn_HighHits * 10.0, drum_flow);
#define DV 0.70710678

// nonlinearity
float nl(float x) {
    return 1.0 / (1.0 + exp(W0 * (W1 * x - 0.5)));
}

vec4 gaussian(vec4 x, vec4 x_nw, vec4 x_n, vec4 x_ne, vec4 x_w, vec4 x_e, vec4 x_sw, vec4 x_s, vec4 x_se) {
    return _G0 * x + _G1 * (x_n + x_e + x_w + x_s) + _G2 * (x_nw + x_sw + x_ne + x_se);
}

bool reset() {
    return reset_sim > 0.5;
}

vec2 normz(vec2 x) {
    return x == vec2(0.0, 0.0) ? vec2(0.0, 0.0) : normalize(x);
}

vec4 advect(vec2 ab, vec2 vUv, vec2 step) {

    vec2 aUv = vUv - ab * ADVECT_DIST * step;

    vec2 n = vec2(0.0, step.y);
    vec2 ne = vec2(step.x, step.y);
    vec2 e = vec2(step.x, 0.0);
    vec2 se = vec2(step.x, -step.y);
    vec2 s = vec2(0.0, -step.y);
    vec2 sw = vec2(-step.x, -step.y);
    vec2 w = vec2(-step.x, 0.0);
    vec2 nw = vec2(-step.x, step.y);

    vec4 u = texture(BuffB, fract(aUv));
    vec4 u_n = texture(BuffB, fract(aUv + n));
    vec4 u_e = texture(BuffB, fract(aUv + e));
    vec4 u_s = texture(BuffB, fract(aUv + s));
    vec4 u_w = texture(BuffB, fract(aUv + w));
    vec4 u_nw = texture(BuffB, fract(aUv + nw));
    vec4 u_sw = texture(BuffB, fract(aUv + sw));
    vec4 u_ne = texture(BuffB, fract(aUv + ne));
    vec4 u_se = texture(BuffB, fract(aUv + se));

    return gaussian(u, u_nw, u_n, u_ne, u_w, u_e, u_sw, u_s, u_se);
}

vec4 advect2(vec2 ab, vec2 vUv, vec2 step, float sc) {

    vec2 aUv = vUv - ab * sc * step;

    // 3x3 neighborhood coordinates
    float step_x = step.x;
    float step_y = step.y;
    vec2 n = vec2(0.0, step_y);
    vec2 ne = vec2(step_x, step_y);
    vec2 e = vec2(step_x, 0.0);
    vec2 se = vec2(step_x, -step_y);
    vec2 s = vec2(0.0, -step_y);
    vec2 sw = vec2(-step_x, -step_y);
    vec2 w = vec2(-step_x, 0.0);
    vec2 nw = vec2(-step_x, step_y);

    vec4 uv = texture(BuffA, fract(aUv));
    vec4 uv_n = texture(BuffA, fract(aUv + n));
    vec4 uv_e = texture(BuffA, fract(aUv + e));
    vec4 uv_s = texture(BuffA, fract(aUv + s));
    vec4 uv_w = texture(BuffA, fract(aUv + w));
    vec4 uv_nw = texture(BuffA, fract(aUv + nw));
    vec4 uv_sw = texture(BuffA, fract(aUv + sw));
    vec4 uv_ne = texture(BuffA, fract(aUv + ne));
    vec4 uv_se = texture(BuffA, fract(aUv + se));

    return _G0 * uv + _G1 * (uv_n + uv_e + uv_w + uv_s) + _G2 * (uv_nw + uv_sw + uv_ne + uv_se);
}

#define SQRT_3_OVER_2 0.86602540378
#define SQRT_3_OVER_2_INV 0.13397459621

vec2 diagH(vec2 x, vec2 x_v, vec2 x_h, vec2 x_d) {
    return 0.5 * ((x + x_v) * SQRT_3_OVER_2_INV + (x_h + x_d) * SQRT_3_OVER_2);
}

vec2 diagV(vec2 x, vec2 x_v, vec2 x_h, vec2 x_d) {
    return 0.5 * ((x + x_h) * SQRT_3_OVER_2_INV + (x_v + x_d) * SQRT_3_OVER_2);
}

vec4 firstPass() {
    vec4 fragColor = vec4(0.0);
    vec2 vUv = _xy.xy / RENDERSIZE.xy;
    vec2 texel = 1. / RENDERSIZE.xy;
    vUv += _uvc*0.0025*stretch;
    vUv += _uvc*0.0025*fisheye*d2;

    if(_mouse.z >= 0.50) {
        vec2 d = _uvc - _muvc * vec2(1.0, 9.0 / 16.0);
        
        vec2 m = 2.0 * normz(d) * exp(-length(d) / 0.03);
        
        vUv += m;
        // div += pow(length(m)*3.0,2.0);
    }

    vec2 n = vec2(0.0, 1.0);
    vec2 ne = vec2(1.0, 1.0);
    vec2 e = vec2(1.0, 0.0);
    vec2 se = vec2(1.0, -1.0);
    vec2 s = vec2(0.0, -1.0);
    vec2 sw = vec2(-1.0, -1.0);
    vec2 w = vec2(-1.0, 0.0);
    vec2 nw = vec2(-1.0, 1.0);

    vec4 u = texture(BuffB, fract(vUv));
    vec4 u_n = texture(BuffB, fract(vUv + texel * n));
    vec4 u_e = texture(BuffB, fract(vUv + texel * e));
    vec4 u_s = texture(BuffB, fract(vUv + texel * s));
    vec4 u_w = texture(BuffB, fract(vUv + texel * w));
    vec4 u_nw = texture(BuffB, fract(vUv + texel * nw));
    vec4 u_sw = texture(BuffB, fract(vUv + texel * sw));
    vec4 u_ne = texture(BuffB, fract(vUv + texel * ne));
    vec4 u_se = texture(BuffB, fract(vUv + texel * se));

    const float vx = 0.5;
    const float vy = SQRT_3_OVER_2;
    const float hx = SQRT_3_OVER_2;
    const float hy = 0.5;

    float di_n = nl(distance(u_n.xy + n, u.xy));
    float di_w = nl(distance(u_w.xy + w, u.xy));
    float di_e = nl(distance(u_e.xy + e, u.xy));
    float di_s = nl(distance(u_s.xy + s, u.xy));

    float di_nne = nl(distance((diagV(u.xy, u_n.xy, u_e.xy, u_ne.xy) + vec2(+vx, +vy)), u.xy));
    float di_ene = nl(distance((diagH(u.xy, u_n.xy, u_e.xy, u_ne.xy) + vec2(+hx, +hy)), u.xy));
    float di_ese = nl(distance((diagH(u.xy, u_s.xy, u_e.xy, u_se.xy) + vec2(+hx, -hy)), u.xy));
    float di_sse = nl(distance((diagV(u.xy, u_s.xy, u_e.xy, u_se.xy) + vec2(+vx, -vy)), u.xy));
    float di_ssw = nl(distance((diagV(u.xy, u_s.xy, u_w.xy, u_sw.xy) + vec2(-vx, -vy)), u.xy));
    float di_wsw = nl(distance((diagH(u.xy, u_s.xy, u_w.xy, u_sw.xy) + vec2(-hx, -hy)), u.xy));
    float di_wnw = nl(distance((diagH(u.xy, u_n.xy, u_w.xy, u_nw.xy) + vec2(-hx, +hy)), u.xy));
    float di_nnw = nl(distance((diagV(u.xy, u_n.xy, u_w.xy, u_nw.xy) + vec2(-vx, +vy)), u.xy));

    vec2 xy_n = u_n.xy + n - u.xy;
    vec2 xy_w = u_w.xy + w - u.xy;
    vec2 xy_e = u_e.xy + e - u.xy;
    vec2 xy_s = u_s.xy + s - u.xy;

    vec2 xy_nne = (diagV(u.xy, u_n.xy, u_e.xy, u_ne.xy) + vec2(+vx, +vy)) - u.xy;
    vec2 xy_ene = (diagH(u.xy, u_n.xy, u_e.xy, u_ne.xy) + vec2(+hx, +hy)) - u.xy;
    vec2 xy_ese = (diagH(u.xy, u_s.xy, u_e.xy, u_se.xy) + vec2(+hx, -hy)) - u.xy;
    vec2 xy_sse = (diagV(u.xy, u_s.xy, u_e.xy, u_se.xy) + vec2(+vx, -vy)) - u.xy;
    vec2 xy_ssw = (diagV(u.xy, u_s.xy, u_w.xy, u_sw.xy) + vec2(-vx, -vy)) - u.xy;
    vec2 xy_wsw = (diagH(u.xy, u_s.xy, u_w.xy, u_sw.xy) + vec2(-hx, -hy)) - u.xy;
    vec2 xy_wnw = (diagH(u.xy, u_n.xy, u_w.xy, u_nw.xy) + vec2(-hx, +hy)) - u.xy;
    vec2 xy_nnw = (diagV(u.xy, u_n.xy, u_w.xy, u_nw.xy) + vec2(-vx, +vy)) - u.xy;

    vec2 ma = di_nne * xy_nne + di_ene * xy_ene + di_ese * xy_ese + di_sse * xy_sse + di_ssw * xy_ssw + di_wsw * xy_wsw + di_wnw * xy_wnw + di_nnw * xy_nnw + di_n * xy_n + di_w * xy_w + di_e * xy_e + di_s * xy_s;

    vec4 u_blur = gaussian(u, u_nw, u_n, u_ne, u_w, u_e, u_sw, u_s, u_se);

    vec4 au = advect(u.xy, vUv, texel);
    vec4 av = advect(u.zw, vUv, texel);

    vec2 dv = av.zw + TIMESTEP * ma;
    vec2 du = au.xy + TIMESTEP * dv;

    if(iMouse.z > 0.0) {
        vec2 d = _xy.xy - iMouse.xy;
        float m = exp(-length(d) / 50.0);
        du += 0.1 * m * normz(d);
    }

    // if (syn_OnBeat > 0.9) {
    //     float m = lum;
    //     du += 0.2 * m ;
    // }

    vec2 init = vec2(_rand(_uv + TIME), _rand(-_uv - TIME));
    // initialize with noise
    if(FRAMECOUNT < 10 || reset()) {
        fragColor = 8.0 * (vec4(-0.5) + vec4(init.xy, init.xy));
    } else {
        vec4 img = _loadUserImage();
        float lum = dot(img.rgb, vec3(1.0)) / 3.0;

        du = length(du) > 1.0 ? normz(du) : du;
        du -= _uvc * pow(_pulse(lum, 0.85, 0.15) * syn_HighHits, 5.0) * media_hits * 3.0;
        du += _uvc * pow(_pulse(lum, 0.85, 0.15) * syn_BassLevel, 5.0) * media_hits * 3.0;
        fragColor = vec4(du, dv);
    }

    return fragColor;
}

vec4 secondPass() {
    vec4 img = _loadUserImage();

    vec4 fragColor = vec4(0.0);
    const float _K0 = -20.0 / 6.0; // center weight
    const float _K1 = 4.0 / 6.0;   // edge-neighbors
    const float _K2 = 1.0 / 6.0;   // vertex-neighbors
    float cs = -3.0 * clouds + 0.5 * snakes;  // curl scale
    const float ls = 3.0;  // laplacian scale
    const float ps = 0.0;  // laplacian of divergence scale
    const float ds = -12.0; // divergence scale
    const float dp = -6.0; // divergence update scale
    const float pl = 0.3;   // divergence smoothing
    const float ad = 6.0;   // advection distance scale
    float pwr = backwards_ripples;  // power when deriving rotation angle from curl
    float amp = 1.0;  // self-amplification
    float upd = smooth_glass;  // update smoothing
    float sq2 = 1.0;  // diagonal weight

    vec2 vUv = _xy.xy / RENDERSIZE.xy;
    vec2 texel = 1. / RENDERSIZE.xy;

    // 3x3 neighborhood coordinates
    float step_x = texel.x;
    float step_y = texel.y;
    vec2 n = vec2(0.0, step_y);
    vec2 ne = vec2(step_x, step_y);
    vec2 e = vec2(step_x, 0.0);
    vec2 se = vec2(step_x, -step_y);
    vec2 s = vec2(0.0, -step_y);
    vec2 sw = vec2(-step_x, -step_y);
    vec2 w = vec2(-step_x, 0.0);
    vec2 nw = vec2(-step_x, step_y);

    vec4 uv = texture(BuffA, fract(vUv));
    vec4 uv_n = texture(BuffA, fract(vUv + n));
    vec4 uv_e = texture(BuffA, fract(vUv + e));
    vec4 uv_s = texture(BuffA, fract(vUv + s));
    vec4 uv_w = texture(BuffA, fract(vUv + w));
    vec4 uv_nw = texture(BuffA, fract(vUv + nw));
    vec4 uv_sw = texture(BuffA, fract(vUv + sw));
    vec4 uv_ne = texture(BuffA, fract(vUv + ne));
    vec4 uv_se = texture(BuffA, fract(vUv + se));

    // uv.x and uv.y are the x and y components, uv.z and uv.w accumulate divergence 

    // laplacian of all components
    vec4 lapl = _K0 * uv + _K1 * (uv_n + uv_e + uv_w + uv_s) + _K2 * (uv_nw + uv_sw + uv_ne + uv_se);

    // calculate curl
    // vectors point clockwise about the center point
    float curl = uv_n.x - uv_s.x - uv_e.y + uv_w.y + sq2 * (uv_nw.x + uv_nw.y + uv_ne.x - uv_ne.y + uv_sw.y - uv_sw.x - uv_se.y - uv_se.x);

    // compute angle of rotation from curl
    float sc = cs * sign(curl) * pow(abs(curl), pwr);

    // calculate divergence
    // vectors point inwards towards the center point
    float div = uv_s.y - uv_n.y - uv_e.x + uv_w.x + sq2 * (uv_nw.x - uv_nw.y - uv_ne.x - uv_ne.y + uv_sw.x + uv_sw.y + uv_se.y - uv_se.x);

    vec2 norm = normz(uv.xy);

    float sdx = uv.z + (dp) * uv.x * div + pl * lapl.z;
    float sdy = uv.w + dp * uv.y * div + pl * lapl.w;

    vec2 ab = advect2(vec2(uv.x, uv.y), vUv, texel, ad).xy;

    // temp values for the update rule
    float ta = amp * ab.x + ls * lapl.x + norm.x * ps * lapl.z + ds * sdx;
    float tb = amp * ab.y + ls * lapl.y + norm.y * ps * lapl.w + ds * sdy;

    // rotate
    float a = ta * cos(sc) - tb * sin(sc);
    float b = ta * sin(sc) + tb * cos(sc);

    vec4 abd = upd * uv + (1.0 - upd) * vec4(a, b, sdx, sdy);

    fragColor = vec4(abd);

    abd.xy = clamp(length(abd.xy) > 1.0 ? normz(abd.xy) : abd.xy, -1.0, 1.0);
    fragColor = vec4(abd);
    return fragColor;

}

// displacement amount
float DISP_SCALE = (-1.0 + 2.0 * switch_direction) * sign(aberration_amount) * pow(aberration_amount, 2.0);

// chromatic dispersion samples
#define SAMPLES 34

// contrast
#define SIGMOID_CONTRAST 8.0

// channels to use for displacement, either xy or zw

vec3 contrast(vec3 x) {
    return 1.0 / (1.0 + exp(-SIGMOID_CONTRAST * (x - 0.5)));
}

/*
    This function supplies a weight vector for each color channel.
    It's analogous to (but not a physically accurate model of)
    the response curves for each of the 3 cone types in the human eye.
    The three functions for red, green, and blue have the same integral
    over [0, 1], which is 1/3.
    Here are some other potential terms for the green weight that 
    integrate to 1/3:
        2.0*(1-x)*x
        10.0*((1-x)*x)^2
        46.667*((1-i)*i)^3
        210.0*((1-x)*x)^4
        924.0*((1-x)*x)^5
    By the way, this series of coefficients is OEIS A004731 divided by 3,
    which is a pretty interesting series: https://oeis.org/A002457
*/
vec3 sampleWeights(float i) {
    return vec3(i * i, 46.6666 * pow((1.0 - i) * i, 3.0), (1.0 - i) * (1.0 - i));
}

vec3 sampleDisp(vec2 uv, vec2 dispNorm, float disp) {
    vec3 col = vec3(0);
    const float SD = 1.0 / float(SAMPLES);
    float wl = 0.0;
    vec3 denom = vec3(0);
    for(int i = 0; i < SAMPLES; i++) {
        vec3 sw = sampleWeights(wl);
        denom += sw;
        if(syn_MediaType > 0.5) {
            col += sw * texture(syn_UserImage, uv + dispNorm * disp * wl).xyz;
        } else {
            vec2 pos = fract(uv * 4.0 * vec2(RENDERSIZE.x / RENDERSIZE.y, 1.0) * vec2(0.145, 2.0) + (dispNorm * disp * wl) * 10.0 * (0.05 + syn_Presence) + vec2(0.0, syn_Time * 0.1));
            pos *= 1.0 - syn_HighHits * 0.5;
            float stepper = 0.1;
            col += sw * vec3((1.0 - step(stepper, pos.x + 0.05) * step(stepper + 0.1, pos.y + 0.05)), (1.0 - step(stepper, pos.x + 0.05) * step(stepper + 0.1, pos.y + 0.05)), (1.0 - step(stepper, pos.x + 0.05) * step(stepper + 0.1, pos.y + 0.05)));
        }
        wl += SD;
    }

    // For a large enough number of samples,
    // the return below is equivalent to 3.0 * col * SD;
    return col / denom;
}

vec4 renderPass() {
    vec4 fragColor = vec4(0.0);
    vec2 texel = 1. / RENDERSIZE.xy;
    vec2 uv = _uv;

    // uv = arbMirror(uv);
    // return vec4(fract(uv.x*10.0),fract(uv.y*10.0), 0.0, 1.0);
    vec2 n = vec2(0.0, texel.y);
    vec2 e = vec2(texel.x, 0.0);
    vec2 s = vec2(0.0, -texel.y);
    vec2 w = vec2(-texel.x, 0.0);

    vec4 dF = texture(BuffA, uv);

    vec4 d_nF = texture(BuffA, fract(uv + n));
    vec4 d_eF = texture(BuffA, fract(uv + e));
    vec4 d_sF = texture(BuffA, fract(uv + s));
    vec4 d_wF = texture(BuffA, fract(uv + w));

    vec2 d = mix(dF.xy, _rotate(dF.xy, dot(dF.zw, vec2(1.0)) * 50.0), distortion + syn_HighLevel * auto_distortion);
    vec2 d_n = mix(d_nF.xy, _rotate(d_nF.xy, dot(d_nF.zw, vec2(1.0)) * 50.0), distortion + syn_HighLevel * auto_distortion);
    vec2 d_e = mix(d_eF.xy, _rotate(d_eF.xy, dot(d_eF.zw, vec2(1.0)) * 50.0), distortion + syn_HighLevel * auto_distortion);
    vec2 d_s = mix(d_sF.xy, _rotate(d_sF.xy, dot(d_sF.zw, vec2(1.0)) * 50.0), distortion + syn_HighLevel * auto_distortion);
    vec2 d_w = mix(d_wF.xy, _rotate(d_wF.xy, dot(d_wF.zw, vec2(1.0)) * 50.0), distortion + syn_HighLevel * auto_distortion); 

    // antialias our vector field by blurring
    vec2 db = 0.4 * d + 0.15 * (d_n + d_e + d_s + d_w);

    float ld = length(db);
    vec2 ln = normz(db);

    if((syn_MediaType > 1.5) && (syn_MediaType < 2.5)) {
        uv.y = 1.0 - uv.y;
    }

    vec3 col = sampleDisp(uv, ln, DISP_SCALE * ld);

    fragColor = vec4(contrast(col), 1.0);
    return fragColor;

}

vec4 renderMain() {
    if(PASSINDEX == 0.0) {
        return firstPass();
    } else if(PASSINDEX == 1.0) {
        return secondPass();
    } else if(PASSINDEX == 2.0) {
        return renderPass();
    }
}
