

			//******** BuffA Code Begins ********
vec2 pixelSize = 1. / RENDERSIZE.xy;
vec2 aspect = vec2(1., RENDERSIZE.y / RENDERSIZE.x);
float a = 1.;
float dist(vec2 p0, vec2 pf) {
    return sqrt((pf.x - p0.x) * (pf.x - p0.x) + (pf.y - p0.y) * (pf.y - p0.y));
}
float d = dist(RENDERSIZE.xy * 0.5, _xy.xy) * (_mouse.x / RENDERSIZE.x + 0.1) * 0.005;
float d2 = dist(RENDERSIZE.xy * 0.5, _xy.xy) * 0.00125;
float o = a * .7071;
vec2 xy_noise = vec2(_noise((_uv + (TIME) * 0.2)), _noise((_uv + (TIME) * 0.2)));
vec2 rot_noise = _rotate(xy_noise, _noise(TIME * 0.1));
vec4 mediaEdges = texture(media_pass_fx, _uv);
vec4 media = texture(media_pass, _uv);
vec4 mixedEdgeCol = mix(media, mediaEdges, edge_mix * 0.95 + 0.1) * media_impact;
vec4 media_influence = mixedEdgeCol * media_impact * 0.25;
float mediahits_mix = mix(media_impact * 0.125, media_impact * 0.125 * (syn_MidLevel * 0.5 * syn_MidPresence + syn_Presence * syn_Level * 0.5), media_hits);
bool mediaOn = _exists(syn_UserImage);
vec4 masked_media = mixedEdgeCol - _invert(_loadMediaAsMask(), 1.0);
float media_lum = PI * 0.25 * sin(length(media_influence * mediahits_mix));
float paintsize = paint_size * 40.;
#define R RENDERSIZE.xy

float sigmoid(float x) {
    return 2. / (1. + exp2(-x)) - 1.;
}
float circle_distance(vec2 uv, vec2 pos, float size, float min) {
    return max(min, 1. - length((uv - pos) * aspect / size));
}

float smooth_circle(vec2 uv, vec2 pos, float size, float ramp) {
    return 0.5 + sigmoid(circle_distance(uv, pos, size, -16.) * ramp) * 0.5;
}

vec2 complex_mul(vec2 factorA, vec2 factorB) {
    return vec2(factorA.x * factorB.x - factorA.y * factorB.y, factorA.x * factorB.y + factorA.y * factorB.x);
}

vec2 complex_div(vec2 numerator, vec2 denominator) {
    return vec2(numerator.x * denominator.x + numerator.y * denominator.y, numerator.y * denominator.x - numerator.x * denominator.y) /
        vec2(denominator.x * denominator.x + denominator.y * denominator.y);
}

vec2 vortex_warp(vec2 uv, vec2 pos, float size, float ramp, vec2 rot) {
    vec2 pos_correct = 0.5 + (pos - 0.5);
    vec2 rot_uv = pos_correct + complex_mul((uv - pos_correct) * aspect, rot) / aspect;
    float smooth_circle = smooth_circle(uv, pos_correct, size, ramp);
    return mix(uv, rot_uv, smooth_circle);
}

vec2 vortex_pair_warp(vec2 uv, vec2 pos, vec2 vel) {
    vec2 aspect = vec2(1., RENDERSIZE.y / RENDERSIZE.x);
    float ramp = 5.;
    float d = 0.2;

    float l = length(vel);
    vec2 p1 = pos;
    vec2 p2 = pos;

    if(l > 0.) {
        vec2 normal = normalize(vel.yx * vec2(-1., 1.)) / aspect;
        p1 = pos - normal * d / 2.;
        p2 = pos + normal * d / 2.;
    }

    float w = l / d * 2.;

  // two overlapping rotations that would annihilate when they were not displaced.

    vec2 circle1 = vortex_warp(uv, p1, d, ramp, vec2(cos(w), sin(w)));
    vec2 circle2 = vortex_warp(uv, p2, d, ramp, vec2(cos(-w), sin(-w)));

    return (circle1 + circle2) / 2.;
}

vec4 D(vec2 U) {
    return texture(BuffD, U / R);
}

vec4 renderPassA() {

    vec4 Q = vec4(0.0);

    vec2 U = _xy;
    vec2 centerModded = _uvc - (Pull - 0.5);

    // U.x*=(1.0-U.y*U.y)*U.x*0.2;
    // U.y=(1.0-U.x*U.x)*U.y*0.2;
    // U = Drift + (U - Drift) * (1.0 - vec2(shock_forward));
    // U = Pull + (U - Pull) * (1.0 - vec2(shock_forward) * aspect * pixelSize * 64.);
    // U += centerModded*length(centerModded)*0.05*(Stretch*Stretch)*sign(Stretch);
    // U += centerModded*(0.5+0.5*length(centerModded))*0.05*shock_back;
    // U += Drift * pixelSize;
    // U -= vec2(Pull.yx - U.yx) * aspect.yx * shear * pixelSize * 128.;
    // U -= vec2(Pull.yx - U.yx) * aspect.yx * pixelSize * 128.;
    vec2 uvc_rot = _rotate(_uvc, TIME);
    U += vec2(rot_noise) * 0.15;
    U += vec2(dot(cos(uvc_rot.x), sin(_uvc.y))) * (Pull + manual_stir*4.) * PI * (2.0 + low) * (1.0 - d2);
    float shock = shock_forward - shock_back * 0.5;
    U -= _uvc * (Zoom) * (1.0 + 2.0 * low + syn_Intensity * 1.5);

    U += _uvc * Stretch * (1.0 + low);
    U += abs(d2) * _uvc * (Fisheye - shock * PI) * PI * (1.0 + 2.0 * low + syn_Intensity * 1.5);

    U += Drift * (1.20 + low);
    U = mix(U, U + vec2(_noise(_uv)) * (warp + manual_stir*2.) * (1. + syn_BassPresence), warp + manual_stir);

    Q = D(U);
    vec2 brushPosition = vec2((sin(TIME * 0.3 + syn_BassTime * 2 * PI * 0.05)), (cos(TIME * 0.25 + syn_BassTime * 2 * PI * 0.04)));

    float twitcher = (pow(syn_BassLevel, 1.3) * 1.25 + manual_stir * 4.);

    U = vortex_pair_warp(U, brushPosition * 0.5 + 0.5, normalize(brushPosition) * (twitcher) * 0.025);

    //float a = TIME+3.*length(U), c = cos(a), s = sin(a);

    float g = smoothTime + 3. * length(U), c = cos(g), s = sin(g);

    mat2 m = mat2(c, -s, s, c);

    vec4 n = D(U + vec2(0, a) * m) + D(U + vec2(a, 0) * m) + D(U + vec2(0, -a) * m) + D(U + vec2(-a, 0) * m) + D(U + vec2(-o, o) * m) + D(U + vec2(o, -o) * m) + D(U + vec2(-o, -o) * m) + D(U + vec2(o, o) * m);

    n *= .125;

    vec4 dx = n - Q;

    Q += dx * vec4(1., .3, 1., 1);

    float x = (.3 + complexity) * Q.x * Q.y * (1. - Q.y);

    Q.y = Q.y + x - 0.025 - .1 * Q.z;

    // Q.x = Q.x + 0.11 * Q.x * (1. - Q.x) - x;
    Q.x = Q.x + (0.1 * (1.0 + mystery * 0.2)) * Q.x * (1. - Q.x) - x;

    Q.z = Q.z * 0.98 + .6 * dx.y + .001 * Q.y + (Split * (1. + syn_BassLevel * 0.1)) * 0.005 + media_lum * 0.1;

    Q = clamp(Q, 0., 1.);

    if(_mouse.z > 0. && length(U - _mouse.xy) < paintsize)
        Q.y = .75;

    // if(_exists(syn_UserImage)) Q.y += _luminance(_edgeDetectSobel(syn_UserImage))*syn_Level*0.2;
    if(_exists(syn_UserImage)) {
        Q.x += media_lum;
        Q.y += media_lum;
    }

    if(FRAMECOUNT <= 1 || Reset != 0.) {

        Q = vec4(1, 0, 0, 0);

        if(length(U - 0.5 * R) < 20.)
            Q.y += 1.;

    }

    return Q;

} 

			//******** BuffB Code Begins ********

#define R RENDERSIZE.xy

vec4 A(vec2 U) {
    return texture(BuffA, U / R);
}

vec4 renderPassB() {

    vec4 Q = vec4(0.0);

    vec2 U = _xy;

    // U += vec2(dot(cos(_uvc.x), sin(_uvc.y))) * Pull * PI*(1.0-d2);

    Q = A(U);

    float g = smoothTime + 3. * length(U), c = cos(g), s = sin(g);

    mat2 m = mat2(c, -s, s, c);

    o = a * .7071;

    vec4 n = A(U + vec2(0, a) * m) + A(U + vec2(a, 0) * m) + A(U + vec2(0, -a) * m) + A(U + vec2(-a, 0) * m) + A(U + vec2(-o, o) * m) + A(U + vec2(o, -o) * m) + A(U + vec2(-o, -o) * m) + A(U + vec2(o, o) * m);

    n *= .125;

    vec4 dx = n - Q;

    Q += dx * vec4(1., .3, 1., 1);

    float x = (.53 + complexity * 0.5) * Q.x * Q.y * (1. - Q.y);

    // if(_exists(syn_UserImage)) {
    //     x = mix(x, sin(_luminance(mixedEdgeCol))*2*PI, syn_MidLevel*0.05*syn_Intensity +syn_Presence*syn_Level*0.05);
    // }

    Q.y = Q.y + x - (0.025 - Thiccness) - .1 * Q.z;
    Q.y *= 1.0 + 0.025 * low;
    Q.z *= 1.0 + low * 0.0125;
    Q.x = Q.x + (0.1 * (1.0 + mystery * 0.2)) * Q.x * (1. - Q.x) - x;
    // Q.x -= media_lum;

    // Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;
    Q.z = Q.z * 0.98 + .6 * dx.y + .001 * Q.y + Split * 0.005 + media_lum * 0.1;

    Q.z *= 1.0 + 0.01 * syn_MidLevel;

    Q = clamp(Q, 0., 1.);

    if(_mouse.z > 0. && length(U - _mouse.xy) < paintsize)
        Q.y = 1. * d2;

    // if(_exists(syn_UserImage)) {
    //     Q.y += media_lum;
    //     Q.x += media_lum;
    // }

    if(FRAMECOUNT <= 1) {

        Q = vec4(1, 0, 0, 0);

        if(length(U - 0.5 * R) < 10.)
            Q.y += 1. * d2;

    }

    return Q;

} 

			//******** BuffC Code Begins ********

#define R RENDERSIZE.xy

vec4 B(vec2 U) {
    return texture(BuffB, U / R);
}

vec4 renderPassC() {

    vec4 Q = vec4(0.0);

    vec2 U = _xy;

    Q = B(U);

    float g = smoothTime + 3. * length(U), c = cos(g), s = sin(g);

    mat2 m = mat2(c, -s, s, c);

    a += 1.0;

    o = a * .7071 + media_lum;

    vec4 n = B(U + vec2(0, a) * m) + B(U + vec2(a, 0) * m) + B(U + vec2(0, -a) * m) + B(U + vec2(-a, 0) * m) + B(U + vec2(-o, o) * m) + B(U + vec2(o, -o) * m) + B(U + vec2(-o, -o) * m) + B(U + vec2(o, o) * m);

    n *= .125;

    vec4 dx = n - Q;

    Q += dx * vec4(1., .3, 1., 1);

    float x = (.3) * Q.x * Q.y * (1. - Q.y);

    Q.y = Q.y + x - (0.025 - Thiccness) - .1 * Q.z;

    // Q.y = Q.y+x-(0.025)-.1*Q.z;

    // Q.x = Q.x + 0.1 * Q.x * (1. - Q.x) - x;
    Q.x = Q.x + (0.1 * (1.0 + mystery * 0.2)) * Q.x * (1. - Q.x) - x;

    Q.z = Q.z * 0.98 + .6 * dx.y + .001 * Q.y + media_lum * 0.1;

    Q = clamp(Q, 0., 2.);

    if(_mouse.z > 0. && length(U - _mouse.xy) < paintsize)
        Q.y = .75;
    if(mediaOn) {
        Q.x += media_lum;
        Q.y += media_lum;
    }

    if(FRAMECOUNT <= 1) {

        Q = vec4(1, 0, 0, 0);

        if(length(U - 0.5 * R) < 10.)
            Q.y += 1.;

    }

    return Q;

} 

			//******** BuffD Code Begins ********

vec4 mediaPass() {
    vec4 media = vec4(0.);
    vec2 U = _xy;
    media = _loadMedia();
    return media;
}
vec4 mediaPassFX() {
    vec4 media_fx = _edgeDetectSobel((media_pass));
    return media_fx;
}
vec4 C(vec2 U) {
    return texture(BuffC, U / R);
}

vec4 renderPassD() {

    vec4 Q = vec4(0.0);

    vec2 U = _xy;

    Q = C(U);

    float g = smoothTime + 3. * length(U), c = cos(g), s = sin(g);

    mat2 m = mat2(c, -s, s, c);

    a += 1.0;

    o = a * .7071 + media_lum;

    vec4 n = C(U + vec2(0, a) * m) + C(U + vec2(a, 0) * m) + C(U + vec2(0, -a) * m) + C(U + vec2(-a, 0) * m) + C(U + vec2(-o, o) * m) + C(U + vec2(o, -o) * m) + C(U + vec2(-o, -o) * m) + C(U + vec2(o, o) * m);

    n *= .125;

    vec4 dx = n - Q;

    Q += dx * vec4(1., .3, 1., 1);

    float x = .3 * Q.x * Q.y * (1. - Q.y);

    x -= Thiccness;
    x -= media_lum * .05;
    //Q.y = Q.y+x-0.025-.1*Q.z;

    // Q.y = Q.y+x-(0.025-Thiccness)-.1*Q.z;
    Q.y = Q.y + x - (0.025) - .1 * Q.z;

    Q.x = Q.x + 0.1 * Q.x * (1. - Q.x) - x;

    Q.z = Q.z * 0.98 + .6 * dx.y + .001 * Q.y + media_lum * 0.1;

    Q = clamp(Q, 0., 1.);

    if(_mouse.z > 0. && length(U - _mouse.xy) < paintsize)
        Q.y = .75;
    if(mediaOn) {
        Q.x += media_lum;
        Q.y += media_lum;
    }

    if(FRAMECOUNT <= 1) {

        Q = vec4(1, 0, 0, 0);

        if(length(U - 0.55 * R) < 10.)
            Q.y += 1. * d2;

    }

    return Q;

}

float ln(vec3 p, vec3 q, vec3 b) {
    return length(p - q - (b - q) * dot(p - q, b - q) / dot(b - q, b - q));
}

vec4 renderMainImage() {

    vec4 Q = vec4(0.0);

    vec2 U = _xy;

    vec4 n = A(U + vec2(0, 1)), e = A(U + vec2(1, 0)), s = A(U - vec2(0, 1)), w = A(U - vec2(1, 0));

    vec2 h = vec2(e.w - w.w, n.w - s.w);

    vec2 g = vec2(e.y - w.y, n.y - s.y);

    vec4 b = A(U);

    vec4 dx = 0.25 * (n + e + s + w) - b;
    vec4 sincol = sin(4. + 1 * (b.z - b.x) * vec4(3, 2, 1, 4));
    // vec4 sincol = sin(4.+.5*(b.z-b.x)*vec4(3,2,1,4));
    // vec4 coscol = cos(6.+.5*(b.z-b.x)*vec4(3,2,1,4));
    vec4 coscol = cos(2. + 1 * (b.x - b.z) * vec4(3, 2, 1, 4));
    vec4 colsmixed = mix(sincol, coscol, dx);
	//Q = (b.w+b.y)*normalize(abs(sin(4.+.5*(b.z-b.x)*vec4(3,2,1,4))));
    Q = (b.w + b.y) * normalize(abs(colsmixed));

    vec3 no = normalize(vec3(g - h, -1));

    vec3 r = reflect(vec3(0, 0, 1), no);

    vec3 l = vec3(2. * R.x, 2. * R.y, 10. * R.y);

    vec3 u = vec3(U, 0);

    vec3 lu = l - u;

    float _o = length(r - lu * dot(r, lu) / dot(lu, lu));
    Q *= 1.0 + media_lum;
    Q = (Q) * (exp(-o) + 2. * exp(-7.5 * o) + 1e2 * exp(-30. * o));

    //Q *= 1.+texture(iChannel1,r);

    Q *= 1. + texture(prism, r.xz);
    Q *= (1.0 + highhits * Flash);

    no = normalize(vec3(g - h, 1));

    r = refract(vec3(0, 0, 1), no, .2);

    l = vec3(3. * R.x, 3. * R.y, 10. * R.y);

    u = vec3(U, 0);

    lu = l - u;

    o = length(r - lu * dot(r, lu) / dot(lu, lu));
    o += media_lum;
    Q = 4. * dx + (Q * 0.9 + 0.1) * (0.4 + exp(-o) + 10. * exp(-9.5 * o));
    Q.rgb = _hueRotate(Q.rgb, hue * 180);
   // Q *= 1.+0.2*sin(Q+3.*texture(iChannel1,r).x*vec4(1,2,3,4));

    //Q *= 1.+0.3*sin(Q+3.*texture(prism,r.yz).x*vec4(1,2,3,4));

    Q *= 1. + 0.4 * sin(Q + 3. * texture(prism, r.yz).x * vec4(1, 2, 3, 4));

   //Q *= 1.+0.2*sin(Q+3.*vec4(1,2,3,4));
    Q += mix(vec4(0.0), masked_media, see_media) * 0.1;

    //Q = mix(Q, masked_media, mediahits_mix);
    Q = mix(Q, _grayscale(_contrast(_brightness(Q, 1.25), 1.5)), BW_Mode);
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

        return mediaPass();

    }
    if(PASSINDEX == 5) {

        return mediaPassFX();

    }
    if(PASSINDEX == 6) {

        return renderMainImage();

    }

}