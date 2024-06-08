//Vial Spiller - Psybernautics - 2023
vec4 b13Battleship(sampler2D image, vec2 uv, vec2 resolution, vec2 direction) {
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
float dist(vec2 p0, vec2 pf) {
    return sqrt((pf.x - p0.x) * (pf.x - p0.x) + (pf.y - p0.y) * (pf.y - p0.y));
}
float cc = dist(RENDERSIZE.xy * 0.5, _xy.xy) * (_mouse.x / RENDERSIZE.x + 0.1) * 0.005;

vec4 si(in vec4 color){
  return vec4(vec3(1.0), 0.0) - color;
}

vec4 sh(float stepx, float stepy, vec2 center, sampler2D tex){
  vec4 tleft = si(texture(tex,clamp(center + vec2(-stepx,stepy), 0.0, 1.0)));
  vec4 left = si(texture(tex,clamp(center + vec2(-stepx,0), 0.0, 1.0)));
  vec4 bleft = si(texture(tex,clamp(center + vec2(-stepx,-stepy), 0.0, 1.0)));
  vec4 top = si(texture(tex,clamp(center + vec2(0,stepy), 0.0, 1.0)));
  vec4 bottom = si(texture(tex,clamp(center + vec2(0,-stepy), 0.0, 1.0)));
  vec4 tright = si(texture(tex,clamp(center + vec2(stepx,stepy), 0.0, 1.0)));
  vec4 right = si(texture(tex,clamp(center + vec2(stepx,0), 0.0, 1.0)));
  vec4 bright = si(texture(tex,clamp(center + vec2(stepx,-stepy), 0.0, 1.0)));
  vec4 x = tleft + 2.0*left + bleft - tright - 2.0*right - bright;
  vec4 y = -tleft - 2.0*top - tright + bleft + 2.0 * bottom + bright;
  vec4 color = sqrt((x*x) + (y*y));
  return color;
}

vec4 eds(sampler2D tex){
    float stepSize = mix(0.2, 1.75, mix(0.0, 0.75+0.125*syn_Level+.125*syn_BassLevel*.05*syn_MidLevel, mix(alt_fb, alt_av, alt_auto)));
    vec2 uv = _uv;
    return sh(stepSize/RENDERSIZE.x, stepSize/RENDERSIZE.y, uv, tex);
}

float ggx3(vec2 uv) {
    return length(texture(h,uv).xyz);
}
    
vec2 ggx2(vec2 uv,float delta) {
    vec2 d=vec2(delta,0);
    return vec2(
        ggx3(uv+d.xy)-ggx3(uv-d.xy),
        ggx3(uv+d.yx)-ggx3(uv-d.yx)
    )/delta;
}

vec2 mc(vec2 uvIn) {
    if (mod(uvIn.x, 2.0) > 1.0){
        uvIn.x = 1.0-uvIn.x;
    }
    if (mod(uvIn.y, 2.0) > 1.0){
        uvIn.y = 1.0-uvIn.y;
    }
    return mod(uvIn, 1.0);
}

vec4 renderMain() {
    if(PASSINDEX == 0){
		return syn_MediaType < 0.5 ? texture(noMedia, _uv) : _loadMedia();
	}
    if(PASSINDEX == 1){
        vec2 uv = _uv;
        uv = mc(uv);
        uv = mix(uv, mix(uv, fract(uv*(1.0+floor(syn_RandomOnBeat*4.5))), .1+floor(syn_ToggleOnBeat*1.5)), tile_media);
        
        float ps = mix(5.0, mix(10.0, 3.0, 1.0-clamp(floor(pow(syn_RandomOnBeat, 2.0)*10.5)*0.1, 0.0, 1.0)), floor(syn_Presence*10.5)*0.1);
        vec2 puv = fract(_pixelate((uv)*ps,0.99)/ps + 0.01*(syn_ToggleOnBeat >= 0.5 ? -1.0 : 1.0));
        uv = mix(uv, puv, mix(pixelate,  syn_BassLevel*syn_BassPresence*0.35+syn_MidLevel*syn_MidPresence*.35+syn_Level*syn_Presence*.25, bass_pix));
        uv = mix(uv, abs(vec2(1.0, 0.0) - uv), mix(0.0, floor(syn_ToggleOnBeat*1.375*syn_Presence), beat_flip));
        vec4 img = texture(a, uv);
		return mix(img, texture(a, uv),tile_media);
	}
	if(PASSINDEX == 2){
        vec2 c_uv = _uvc;
        c_uv -= vec2(sin(bass_time), cos(TIME));
        c_uv = _rotate(c_uv, bass_time);
        
        float s = mix(0.0025, 0.125, _size);
        float squareBlast = _sqPulse(abs(c_uv.x + _fbm(_uv*distort*.5+ length(_uv)*0.1*distort+mid_time*0.1)), 1.0-syn_MidLevel, s*mix(0.2, 1.0, _size)+s*syn_MidLevel);
        float p_val = squareBlast ;
        vec3 noise_text = _hsv2rgb(vec3(_noise(p_val))+vec3(p_val+bass_count, 0.0, 0.0));
        vec4 col =  vec4(vec3(noise_text.x), 1.0);
        vec4 lolText = texture(a, _uv);
        // return mix(vec4(0.0), abs(texture(b, _uv)-texture(bp, _uv)), keep_color);
        if (_mouse.z != 0.){

        }
		return  disruptor < 0.5 ? mix(vec4(0.0), mix(texture(b, _uv), abs(texture(b, _uv)-texture(bp, _uv)), mix(alt_fb, alt_av, alt_auto)), keep_color) : col*mix(lolText, lolText+abs(texture(b, _uv)), keep_color);
	}
    else if (PASSINDEX == 3) {
        vec4 temp = vec4(0.0);
        vec4 img = texture(c,_uv);
        vec3 imgHSV = _rgb2hsv(img.xyz);
        vec2 viewPort = (vec2(0.0, 0.0))*mix(-0.005, 0.005, imgHSV.x);
        vec2 uwu = ( ( _uv - 0.5 ) / ( 1.0 - mix(0.00, 0.000625, syn_BassPresence*0.75+syn_MidPresence*0.25)*(1.0 - show_media)*syn_Presence) + 0.5 + viewPort);
        vec2 uwuSave = ( ( _uv - 0.5 ) / ( 1.0 + mix(0.005, -0.015, imgHSV.x)*0.0) + 0.5 + viewPort);
        _uv2uvc(uwu);
        uwu += _uvc*0.01*zoom;
        uwu += _uvc*0.04*fisheye*cc;
        vec2 uwuFluid = uwu;
        vec2 uwuFluidSave = uwu;
        float size = mix(0.005, mix(0.0125, 0.025, mix(alt_fb, alt_av, alt_auto)), syn_Presence);
        float flow = sin(uwu.y+TIME);
        uwu.x += img.r*size + flow;
        uwu.y += img.g*size + flow;
        uwu -= img.b*size + flow;
        float offset = 1.0/RENDERSIZE.x;
        vec4 p = texture(f, (uwuSave + vec2(0.0, 0.0)*offset));
        vec4 n = texture(bp, (uwuSave + vec2(0.0, 1.0)*offset));;
        vec4 e = texture(bp, (uwuSave + vec2(1.0, 0.0)*offset));;
        vec4 s = texture(bp, (uwuSave + vec2(0.0, -1.0)*offset));;
        vec4 w = texture(bp, (uwuSave + vec2(-1.0, 0.0)*offset));;
        vec4 ne = texture(bp, (uwuSave + vec2(1.0, 1.0)*offset));;
        vec4 nw = texture(bp, (uwuSave + vec2(-1.0, 1.0)*offset));;
        vec4 se = texture(bp, (uwuSave + vec2(1.0, -1.0)*offset));;
        vec4 sw = texture(bp, (uwuSave + vec2(-1.0, -1.0)*offset));;
        vec4 media = texture(b, _uv);
        media.xyz = media.xyz;

        vec4 fluidFB = (n+e+s+w+ne+nw+se+sw);
        vec3 selfEdge = _rgb2hsv(p.xyz);
        selfEdge = mix(selfEdge, (_rgb2hsv(media.xyz)+selfEdge)/2.0, clamp(1.0, 0.0,1.0));
        vec3 fluidFBHSV = _rgb2hsv(fluidFB.xyz/mix(12.0, mix(4.0, 4.0, clamp(1.0, 0.0,1.0)), selfEdge.z));
        fluidFB /= mix(12.0, 2.0, clamp(1.0 - selfEdge.z, 0.0, 1.0));
        float fluidSize = mix(size, offset,fluidFBHSV.z);
        float fs = fluidSize*mix(offset, selfEdge.z, selfEdge.z);
        uwuFluid.x -= fluidFB.r*fluidSize + syn_MidHighPresence*fs;
        uwuFluid.y -= fluidFB.g*fluidSize + syn_BassPresence*fs;
        uwuFluid += fluidFB.b*fluidSize + syn_BassHits*fs;
        uwu = mix(uwuFluid, uwu, 0.0);
        vec4 fragColor = vec4(0.0);
        _uvc2uv(uwu);
        vec2 uwu2 = uwu;
        vec2 uwu3 = uwu;
        float gp = _glitch*(syn_Presence*(syn_BassLevel*0.01 + syn_HighLevel*0.01 + syn_MidHighLevel*0.01 + + syn_MidLevel*0.01)*0.5);
        float glitchSize = pow(2.0, pow(2.0, mix(1.0, 2.5, (floor(syn_RandomOnBeat*10.5)*0.1)*_glitch)))*0.5;
        float shift = fract(uwu.y*glitchSize)*gp;
        float shift_switch = mix(-1.0, 1.0, floor(syn_ToggleOnBeat*1.5));
        shift *= shift_switch;
        float fgs = floor(glitchSize);
        float shiftMod_x = mix(sin(uwu2.x*fgs + gp), 2.0*glitch_style, glitch_style);
        float shiftMod_y = mix(sin(uwu2.y*fgs + gp), 2.0*glitch_style, glitch_style);
        if (mod(fract(uwu.y*glitchSize), 2.0) > 0.0) {
            uwu2 += vec2(shift*shiftMod_y, mix(0.0, 0.5, floor(syn_RandomOnBeat*1.5)));
        }
        else {
            uwu2 -= vec2(shift*shiftMod_y, mix(0.0, 0.5, floor(syn_RandomOnBeat*1.5)));
        }
        uwu = mix(uwu, uwu2, clamp(_glitch, 0.0, 1.0));
        shift = fract(uwu.x*glitchSize)*gp;
        if (mod(fract(uwu.x*glitchSize), 2.0) >0.0) {
            uwu3 += vec2(0.0, shift*shiftMod_x);
        }
        else {
            uwu3 -= vec2(0.0, shift*shiftMod_x);
        }
        uwu = mix(uwu, uwu3, clamp(_glitch,0.0,1.0));
        uwu = mix(uwu2, uwu3, floor(syn_RandomOnBeat+0.5));
        float thresh = mix(0.0,mix(kp_thresh, mix(0.25, 0.9, syn_BassLevel), bass_thresh), keep_color);
        vec3 ihsv = _rgb2hsv(img.xyz);   
        if (img.x <= thresh && img.y <= thresh && img.z <= thresh) {
            vec4 fb = texture(ep, uwu);
            fb = mix(fb, abs(fb+media)/2.0,mix(show_media, show_media , mix(alt_fb, alt_av, alt_auto)));
            fb.xyz *= mix(1.0, 1.001, 1.0 - syn_BassPresence);
            fb.xyz = _rgb2hsv(fb.xyz);
            fb.y += mix(0.0, 0.0005, fb.x);
            fb.z = clamp(fb.z, 0.0, 1.0);
            fb.xyz = _hsv2rgb(fb.xyz); 
            img = mix(img, fb, traceMix);
        }
        else {
            img = texture(a, _uv);
        }
        return clamp(img, 0.0, 1.0);
    }
    else if(PASSINDEX == 4.0){
        return texture(d, _uv);
    }
    else if(PASSINDEX == 5.0){ 
        vec4 e = eds(d);
        vec4 m = texture(b, _uv);
        return mix(abs(m-e), e*1.0, mix(alt_fb, alt_av, alt_auto));
    }
    else if (PASSINDEX == 6.0) {
        return b13Battleship(f, _uv, RENDERSIZE, _uvc);
    }
    else if(PASSINDEX == 7.0){ 
        vec4 ogColor = texture(ep, _uv);
        vec2 uwu = ( ( _uv - 0.5 ) / ( 1.0 + mix(0.0, 0.05, 0.0)) + 0.5);
        vec2 uv = _uv;
        vec4 feedback = texture(ep, uv);
        vec4 edge = texture(h, uv);
        vec4 fragColor = vec4(1.0);
        vec3 e = _rgb2hsv(edge.xyz);
        vec4 blur = texture(bp, _uv);
        float thresh = 0.05;
        if (ogColor.x <= thresh && ogColor.y <= thresh && ogColor.z <= thresh) {
            feedback = mix(feedback, clamp(feedback+blur*3.0, 0.0, 1.0), 1.0-e.z);
        }
        feedback = clamp(feedback, 0.0, 1.0);
        return mix(ogColor, feedback, traceMix);
    }
    if(PASSINDEX == 8.0){
        vec4 img = mix(texture(g, _uv), abs(texture(g, _uv) - eds(g)*mix(1.0, 2.0, mix(alt_fb, alt_av, alt_auto))), 1.0);
        return img;
	}
    
    if(PASSINDEX == 9.0) {
        vec2 uv = _uv;
        vec3 n = vec3(ggx2(uv,1.0/RENDERSIZE.y), 50.0);
        n *= n;
        n=normalize(n);
        vec4 r=vec4(n,1);
        vec3 lt = normalize(vec3(1.0,1.0,1.5));
        float df=clamp(dot(n,lt),0.0,1.0);
        float sp=clamp(dot(reflect(lt,n),vec3(0,0,mix(-1, -20, csf))),0.0,1.0);
        sp=pow(sp,mix(5.0, 2.0, csf))*mix(1.5, 0.5, csf);
        vec4 lts = vec4(df)+sp;
        vec4 pp = texture(h, _uv);
        r = pp*mix(vec4(1.0), lts, textured_feedback);
        return r;
    }
    if(PASSINDEX == 10.0) {
       vec4 iSave = texture(jp, _uv);
        vec4 img = iSave*mix(2.5, 3.25, mix(am_i_tenorless, syn_Presence+am_i_tenorless, xxxtra));
        for (float i = 0.0; i < 4.0; i++) {
            img = abs(img - vec4(1.0));
        }
        return _mix3(iSave, (iSave*img), img, mix(am_i_tenorless, syn_Presence+am_i_tenorless, xxxtra));
        return texture(jp, _uv);
    }
}