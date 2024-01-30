//G3LT4B2 - Fine Control - Psybernautics (Alex Tiemann) - 2023

float smin(float a, float b, float k) {
    float h = max(k - abs(a - b), 0.0);
    return min(a, b) - h * h * 0.25 / k;
}

vec2 rotateCenter(vec2 uvIn, float amount) {
    uvIn.y += (RENDERSIZE.x - RENDERSIZE.y) / RENDERSIZE.x;
    uvIn *= vec2(1.0, RENDERSIZE.y / RENDERSIZE.x);
    _uv2uvc(uvIn);
    uvIn = _rotate(uvIn, amount);
    _uvc2uv(uvIn);
    uvIn /= vec2(1.0, RENDERSIZE.y / RENDERSIZE.x);
    uvIn.y -= (RENDERSIZE.x - RENDERSIZE.y) / RENDERSIZE.x;
    return uvIn;
}

vec3 _grad3(vec3 col1, vec3 col2, vec3 col3, float mixVal) {
    mixVal *= 2.0;
    float mix1 = clamp(mixVal, 0.0, 1.0);
    float mix2 = clamp(mixVal - 1.0, 0.0, 1.0);
    return mix(mix(col1, col2, mix1), mix(col2, col3, mix2), step(1.0, mixVal));
}

vec4 texMirror(sampler2D samplerIn, vec2 uvIn) {

    if(mod(uvIn.x, 2.0) > 1.0) {
        uvIn.x = 1.0 - uvIn.x;
    }
    if(mod(uvIn.y, 2.0) > 1.0) {
        uvIn.y = 1.0 - uvIn.y;
    }
    return texture(samplerIn, mod(uvIn, 1.0));

}

vec4 blur13(sampler2D image, vec2 uv, vec2 resolution, vec2 direction) {

    vec4 color = vec4(0.0);

    vec2 off1 = vec2(1.411764705882353) * direction;

    vec2 off2 = vec2(3.2941176470588234) * direction;

    vec2 off3 = vec2(5.176470588235294) * direction;

    color += texture(image, uv) * 0.1964825501511404;

    color += texture(image, uv + (off1 / resolution)) * 0.2969069646728344;

    color += texture(image, uv - (off1 / resolution)) * 0.2969069646728344;

    color += texture(image, uv + (off2 / resolution)) * 0.09447039785044732;

    color += texture(image, uv - (off2 / resolution)) * 0.09447039785044732;

    color += texture(image, uv + (off3 / resolution)) * 0.010381362401148057;

    color += texture(image, uv - (off3 / resolution)) * 0.010381362401148057;

    return color;

}

vec4 sobelIntensity(in vec4 color) {

    return color * 0.5;

}

vec4 sobelHelper(float stepx, float stepy, vec2 center, sampler2D tex) {

    vec4 tleft = sobelIntensity(texture(tex, clamp(center + vec2(-stepx, stepy), 0.0, 1.0)));

    vec4 left = sobelIntensity(texture(tex, clamp(center + vec2(-stepx, 0), 0.0, 1.0)));

    vec4 bleft = sobelIntensity(texture(tex, clamp(center + vec2(-stepx, -stepy), 0.0, 1.0)));

    vec4 top = sobelIntensity(texture(tex, clamp(center + vec2(0, stepy), 0.0, 1.0)));

    vec4 bottom = sobelIntensity(texture(tex, clamp(center + vec2(0, -stepy), 0.0, 1.0)));

    vec4 tright = sobelIntensity(texture(tex, clamp(center + vec2(stepx, stepy), 0.0, 1.0)));

    vec4 right = sobelIntensity(texture(tex, clamp(center + vec2(stepx, 0), 0.0, 1.0)));

    vec4 bright = sobelIntensity(texture(tex, clamp(center + vec2(stepx, -stepy), 0.0, 1.0)));

    vec4 x = tleft + 2.0 * left + bleft - tright - 2.0 * right - bright;

    vec4 y = -tleft - 2.0 * top - tright + bleft + 2.0 * bottom + bright;

    vec4 color = sqrt((x * x) + (y * y));

    return color;

}

vec4 edgeDetectSobel(sampler2D tex) {

    float z = _rgb2hsv(texture(glitchPass, _uv).xyz).x;

    float stepSize = mix(0.0, 1.0, mix(1.0 - z, z, 0.0));

    vec2 uv = _uv;

    return sobelHelper(stepSize / RENDERSIZE.x, stepSize / RENDERSIZE.y, uv, tex);

}

float getVal(vec2 uv) {
    return length(texture(finalMirrPass, uv).xyz);
}

vec2 getGrad(vec2 uv, float delta) {
    vec2 d = vec2(delta, 0);
    return vec2(getVal(uv + d.xy) - getVal(uv - d.xy), getVal(uv + d.yx) - getVal(uv - d.yx)) / delta;
}

vec2 mirrorCoords(vec2 uvIn) {

    if(mod(uvIn.x, 2.0) > 1.0) {
        uvIn.x = 1.0 - uvIn.x;
    }

    if(mod(uvIn.y, 2.0) > 1.0) {
        uvIn.y = 1.0 - uvIn.y;
    }

    return mod(uvIn, 1.0);
}

vec3 getColorGrad(vec3 a, vec3 b, float f) {
    return _grad3(a, b, b, fract((length(a) - bass_time * 0.1) * f));
}

vec4 getFBSum(sampler2D image, vec2 uv, float offset) {

    vec4 n = texture(image, (uv + vec2(0.0, 1.0) * offset));
    ;

    vec4 e = texture(image, (uv + vec2(1.0, 0.0) * offset));
    ;

    vec4 s = texture(image, (uv + vec2(0.0, -1.0) * offset));
    ;

    vec4 w = texture(image, (uv + vec2(-1.0, 0.0) * offset));
    ;

    vec4 ne = texture(image, (uv + vec2(1.0, 1.0) * offset));
    ;

    vec4 nw = texture(image, (uv + vec2(-1.0, 1.0) * offset));
    ;

    vec4 se = texture(image, (uv + vec2(1.0, -1.0) * offset));
    ;

    vec4 sw = texture(image, (uv + vec2(-1.0, -1.0) * offset));
    ;

    return (n + e + s + w + ne + nw + se + sw);
}

vec4 fisheyeEffect(sampler2D text, vec2 uv, float strength) {
    // Transform the uv coordinates to the [-1, 1] range
    vec2 centeredUV = (uv * 2.0) - vec2(1.0);

    // centeredUV = ( ( centeredUV - 0.5 ) / ( 1.0 ) + 0.5 );

    // Calculate the distance from the center
    float dist = length(centeredUV);

    // Apply the fisheye distortion
    float distortion = pow(dist, strength);

    // float zoomFactor = 1.0 + (strength);
    //  distortion = 1.0 / (dist * zoomFactor + 1.0);

    // Transform the distorted coordinates back to the [0, 1] range
    vec2 distortedUV = (centeredUV * distortion + vec2(1.0)) * 0.5;

    // Sample the texture with the distorted coordinates
    return texture(text, distortedUV + vec2(-TIME * 0.0, 0.0));
}

vec4 renderMain() {
    if(PASSINDEX == 0) {
        return _loadUserImage();
    }

    if(PASSINDEX == 1) {
        vec2 uv = _uv;
        uv = mirrorCoords(uv);
        uv = mix(uv, mix(uv, fract(uv * (1.0 + floor(syn_RandomOnBeat * 2.5))), floor(syn_RandomOnBeat * 10.5) * 0.1), tile_media);
        vec4 m = texture(mediaPass, rotateCenter(uv, mix(0.0, mix(0.0, mix(-0.01, -0.1, syn_BassPresence), abs(sin(uv.x * cos(uv.y * mix(-twist, twist, syn_BassPresence)) * mix(twist, -twist, syn_BassPresence) * length(_uvc)))), twist * 0.2) * (0.1 - 0.01 * twist))) * vec4(vec3(1.5), 1.0);
        return m;
    }

    if(PASSINDEX == 2) {
        //mix(texture(mediaHelper, _uv), blur13(mediaHelper, uwu, RENDERSIZE, vec2(0.5)), clamp(1.0-0.9975+0.5, 0.0, 1.0))
        vec4 f = FRAMECOUNT < 3 ? texture(noise, _uv) : abs(edgeDetectSobel(rVertMirrPass)) * 0.95;//
        // f.a *= 0.05;
        return f;
    } else if(PASSINDEX == 3) {
        vec4 temp = texture(postFXPass, _uv);
        vec4 img = mix(temp, texture(mediaPassFX, _uv), 0.0);// show_media);

        if(traceMix < 0.1) {
            return img;
        }

        float border = abs(_uvc.x) >= 1.0 ? -1.0 : 1.0;

        vec2 viewPort = vec2(syn_BassPresence, syn_MidPresence) * mix(-0.001, 0.001, 0.0);
        vec2 uwu = ((_uv - 0.5) / (1.0 + mix(0.000, mix(0.0, 0.0025 + syn_BassPresence * 0.005, syn_Presence), zoom) * border) + 0.5 - _uvc * 0.0025 * syn_BassPresence * border);
        vec2 uwuSave = ((_uv - 0.5) / (1.0 + 0.00) + 0.5);

        // img.xyz = texture(finalMirrPass, uwu).xyz;

        _uv2uvc(uwu);
        vec2 uwuFluid = uwu;
        vec2 uwuFluidSave = uwu;
        // float size = 0.01;
        // float flow = sin(uwu.y);

        // uwu.x += img.r*size + flow;

        // uwu.x -= img.g*size + flow;

        // uwu -= img.b*size + flow;

        float offset = 20.0 / RENDERSIZE.x;

        vec4 p = texture(secondPassB, (_uv + vec2(0.0, 0.0) * offset));

        vec4 media = texture(mediaPassFX, _uv);

        // media.xyz = media.xyz;

        // vec3 bwMedia = vec3(dot(media.xyz,media.xyz));

        vec4 fluidFB = getFBSum(secondPassB, uwuSave, offset);

        vec3 selfEdge = _rgb2hsv(p.xyz);

        // selfEdge = mix(selfEdge, (_rgb2hsv(media.xyz)+selfEdge)/2.0, clamp(show_media+hint_media*0.21, 0.0,1.0));

        // vec3 fluidFBHSV = _rgb2hsv(fluidFB.xyz/mix(8.0, mix(12.0, 12.0, clamp(show_media+hint_media*0.21, 0.0,1.0)), selfEdge.z));

        fluidFB /= mix(12.0, 4.0, mix(flow, syn_BassLevel, auto_flow));// mix(7.0, mix(6.0, 6.0, clamp(show_media+hint_media*0.21, 0.0,1.0)), clamp(mix(fluidFBHSV.z, 1.0 - fluidFBHSV.z, selfEdge.z), 0.0, 1.0));

        float fluidSize = mix(0.125, offset, selfEdge.z) * border;
        // float fluidSize = mix(mix(flow, mix(-0.25, 0.25, syn_RandomOnBeat), auto_flow), offset,selfEdge.z);

        // if (abs(_uvc.y) >= 0.9 || abs(_uvc.x) >= 0.9){
        //     fluidSize *= -1.0;
        // }

        uwuFluid.y += fluidFB.r * fluidSize;// - syn_MidHighPresence*fs;

        // - syn_BassPresence*fs;

        uwuFluid.y -= fluidFB.b * fluidSize;// + syn_HighPresence*fs;
        uwuFluid.x -= fluidFB.g * fluidSize * sin(_uvc.y * 2.5);
        uwu = mix(uwuFluid, uwuFluidSave, 0.0);

        vec4 fragColor = vec4(0.0);

        _uvc2uv(uwu);

        float thresh = 0.95;//clamp(show_media*2.0, 0.0, 0.5);
        vec3 ihsv = _rgb2hsv(img.xyz);
        if(img.x <= thresh || img.y <= thresh || img.z <= thresh) {
            //uwu.x = mix(uwu.x, 1.0-uwu.x, mix(0.0, fb_fckry, fb_fckry));
            vec4 fb = texture(secondPassB, rotateCenter(uwu, mix(0.0, mix(tS, -tS, tan(sin(uwu.x * cos(uwu.y * mix(-5.0, 5.0, syn_BassPresence)) * mix(5.0, -5.0, syn_BassPresence) * length(uwu)))), twist) * (0.1 - 0.05 * twist)));
            // fb = mix(fb, mix(fb,media,0.5),show_media);
            // fb.xyz *= mix(1.0, 1.001, 1.0 - syn_BassPresence);
            vec3 mediaHSV = _rgb2hsv(media.xyz);
            fb.xyz = _rgb2hsv(fb.xyz);
            fb.x += mix(0.0, 0.05, fb.z);
            // fb.x = mix(fb.x, mediaHSV.x,show_media);
            // fb.z = clamp(fb.z, 0.0, 0.9);
            fb.xyz = _hsv2rgb(fb.xyz);
            img = mix(img, fb, traceMix);
        }

        return clamp(img, 0.0, 1.0);
    } else if(PASSINDEX == 4.0) {

        return texture(secondPass, _uv);
    } else if(PASSINDEX == 5.0) {

        vec4 es = edgeDetectSobel(secondPass);
        vec4 e = mix(es, abs(edgeDetectSobel(mediaPassFX) - es * edge_only), show_media);

        vec4 m = mix(texture(secondPass, _uv), texture(mediaPassFX, _uv), 0.0);

        return abs(m - e);
    } else if(PASSINDEX == 6.0) {

        vec2 uwu = ((_uv - 0.5) / (1.0 - mix(0.001, 0.05, 0.0)) + 0.5 - _uvc * 0.001);

        float offset = 2.0 / RENDERSIZE.x;
        vec4 fluidFB = getFBSum(rVertMirrPass, uwu, offset);
        fluidFB /= mix(12.0, 2.0, mix(flow, pow(syn_Presence * syn_Intensity, 2.0), auto_flow));

        uwu.x += fluidFB.r * offset;
        uwu.x -= fluidFB.b * offset;
        uwu.y -= fluidFB.g * offset;
        uwu.y += fluidFB.a * offset;

        vec4 ogColor = texture(postFXPass, uwu);

        vec4 feedback = texture(secondPass, _uv);
        //fragColor.xyz = mix(fragColor.xyz, reflect(fragColor.xyz, normalize(vec3(_uv, 10.0) )), 0.25);
        vec4 edge = texture(finalMirrPass, _uv);
        vec4 fragColor = vec4(1.0);
        vec3 e = _rgb2hsv(edge.xyz);
        vec2 uwu2 = uwu;
        _uv2uvc(uwu2);
        vec4 blur = blur13(edgePass, _uv, RENDERSIZE, uwu2 * 5.0);
        float thresh = 0.95;

        if(feedback.x <= thresh || feedback.y <= thresh || feedback.z <= thresh) {

            float swearl = mix(-0.0025, 0.0025, mix(swirl, syn_Level * syn_ToggleOnBeat, auto_swirl));
            feedback = texture(rVertMirrPass, rotateCenter(uwu, sin(uwu2.y) * cos(uwu2.x) * swearl - swearl));
            feedback = mix(feedback, abs(texture(mediaPassFX, uwu) - feedback * diff_media), show_media * edge_only);
            feedback = mix(feedback, clamp(feedback + blur * 3.0, 0.0, 1.0), 1.0 - e.z);
            feedback.xyz = _rgb2hsv(feedback.xyz);
            // feedback.z -= 0.005;
            feedback.x += 0.001;
            feedback.x = clamp(feedback.x, 0.0, 1.0 - palette);

            feedback.xyz = _hsv2rgb(feedback.xyz);
        }

        feedback = clamp(feedback, 0.0, 1.0);
        return mix(ogColor, feedback, traceMix);
    }

    if(PASSINDEX == 7.0) {
        vec4 img = blur13(glitchPass, _uv, RENDERSIZE, vec2(_uvc * -1.0));

        img = fisheyeEffect(glitchPass, _uv, distort);
        return img;
    }

    if(PASSINDEX == 8.0) {

        vec2 uv = _uv;
        vec2 uv2 = _uv;

        uv2.x = 1.0 - uv2.x;

        vec4 i = texture(rVertMirrPass, uv);
        vec4 i2 = texture(rVertMirrPass, uv2);

        i = mix(i, abs(i + i2) * 0.5, tunnel);

        return i;
    }

    if(PASSINDEX == 9.0) {

        vec2 uv = _uv;
        vec3 n = vec3(getGrad(uv, 1.0 / RENDERSIZE.y), 20.25);

        n *= n;
        n = normalize(n);
        vec4 regColor = vec4(n, 1);
        vec3 light = normalize(vec3(1.0, 1.0, 1.5));
        float diff = clamp(dot(n, light), 0.5, 1.0);
        float spec = clamp(dot(reflect(light, n), vec3(0, 0, -1)), 0.0, 1.0);
        spec = pow(spec, mix(64.0, 2.0, 0.0)) * mix(.5, 1.0, 0.0);
        vec4 lights = vec4(diff) + spec;
        regColor = texture(finalMirrPass, _uv);

        regColor = abs(abs(regColor * 2.0 - vec4(1.0)) - vec4(1.0));

        regColor *= mix(lights, vec4(1.0), 0.0);
        return regColor;
    }
}