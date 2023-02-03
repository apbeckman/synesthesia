

//CC0 1.0 Universal https://creativecommons.org/publicdomain/zero/1.0/
//To the extent possible under law, Blackle Mori has waived all copyright and related or neighboring rights to this work.

const float bpm = 168.;

vec3 erot(vec3 p, vec3 ax, float ro) {
    return mix(dot(p,ax)*ax,p,cos(ro))+sin(ro)*cross(ax,p);
}

vec3 arot(vec3 p, vec3 a, vec3 b) {
    vec3 ax = normalize(cross(a,b));
    float ro = acos(dot(a,b));
    return erot(p,ax,ro);
}

vec3 urot(vec3 p, vec3 a, vec3 b, vec3 u) {
    float ang = atan(dot(u,cross(a,b)), dot(a,b));
    p = erot(p, u, ang);
    a = erot(a, u, ang);
    return arot(p,a,b);
}

vec3 face(vec3 p) {
    vec3 ap = abs(p); vec4 k = vec4(sign(p),0);
    if (ap.x >= max(ap.y,ap.z)) return k.xww;
    if (ap.y >= max(ap.x,ap.z)) return k.wyw;
    if (ap.z >= max(ap.y,ap.x)) return k.wwz;
}

vec3 edge(vec3 p) {
    vec3 mask = vec3(1)-abs(face(p));
    vec3 v = sign(p);
    vec3 a = v*mask.zxy;
    vec3 b = v*mask.yzx;
    return distance(p,a)<distance(p,b)?a:b;
}


#define FK(k) floatBitsToInt(cos(k))^floatBitsToInt(k)
float hash(float a, float b) {
    int x = FK(a),y=FK(b);
    return float((x*x-y)*(y*y+x)+x)/2.14e9;
}

float box(vec2 p, vec2 d) {
    vec2 q = abs(p)-d;
    return length(max(q,0.)) + min(0.,max(q.x,q.y));
}

float linedist(vec2 p, vec2 a, vec2 b) {
    float k = dot(p-a,b-a)/dot(b-a,b-a);
    return distance(p,mix(a,b,clamp(k,0.,1.)));
}

float pentagram(vec2 p) {
    float dist = 10000.;
    for (int i = 0; i < 5; i ++ ){
        float ang = float(i)*3.14/5.*2.;
        float ang2 = float(i+2)*3.14/5.*2.;
        vec2 a = vec2(sin(ang), cos(ang));
        vec2 b = vec2(sin(ang2), cos(ang2));
        dist = min(linedist(p,a,b),dist);
    }
    return dist;
}

vec3 lpos;
vec3 id;
float pedistal;
float l1;
float dbg;
float torch;
float scene(vec3 p) {
    vec3 op = p;
    
    vec2 pc = vec2(length(p.xy),p.z);
    pedistal = box(pc+vec2(0.,3.8), vec2(2.,0.1));
    float ang = atan(p.x,p.y);
    pedistal += smoothstep(-.2,.2,sin(ang*100.))*.02;
    pedistal = min(pedistal, box(pc+vec2(0.,3.8), vec2(1.9,0.15)));
    float ang2 = round(ang*4.)/4.;
    vec2 clos = vec2(sin(ang2),cos(ang2))*1.7;
    
    float ang3 = round(ang*2.)/2.;
    vec2 clos2 = vec2(sin(ang3),cos(ang3))*2.;
    
    pc = vec2(length(p.xy-clos),p.z);
    vec2 pc2 = vec2(length(p.xy-clos2),p.z);
    dbg = sin(min(pentagram(p.xy),pc.x)*100.);
    pedistal = min(pedistal, box(pc+vec2(0.,4.), vec2(.01,0.8)))-0.03;
    pedistal = min(pedistal, box(pc2+vec2(0.,6.5), vec2(.1,2.8)));
    torch = box(pc+vec2(0.,3.2), vec2(.05,0.05));
    pedistal = min(torch, pedistal);
    
    
    vec3 beyp = vec3(7,0,0);
    id = floor(p)+.5;
    vec3 m = sign(mod(id,2.)-1.);
    if (m.x*m.y*m.z < 0.) id += face(p-id);
    float sd = hash(hash(id.x,id.z),id.y);
    if (length(id) < 6. || sd < -0.5 || distance(id,beyp) < 3.7) {
        id += edge(p-id);
    }
    p -= id;
    float balls = max(5.-length(op),length(p)-.7);
    l1 = distance(op,lpos) - .3;
    if (distance(op,beyp) < 3.1) {
        id = beyp;
    }
    balls = min(distance(op,beyp) - 3.,balls);
    return min(pedistal, min(l1, balls));
}

vec3 norm(vec3 p) {
    mat3 k = mat3(p,p,p)-mat3(0.01);
    return normalize(scene(p)-vec3(scene(k[0]),scene(k[1]),scene(k[2])));
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-.5*RENDERSIZE.xy)/RENDERSIZE.y;

    float bpmt = TIME/60.*bpm;
    float t = mix(pow(sin(fract(bpmt)*3.14/2.),50.) + floor(bpmt), bpmt, 0.8);
    lpos = vec3(sin(t),cos(t),sin(t/3.))*4.;
    vec3 cam = normalize(vec3(.8 + sin(t*3.14)*.2,uv));
    vec3 init = vec3(-3.5+cos(t*3.14/2.),0,0) + cam*.1;
    
    cam = erot(cam, vec3(0,1,0), sin(t/2.-.5)*.5);
    init = erot(init, vec3(0,1,0), sin(t/2.)*.5);
    cam = erot(cam, vec3(0,0,1), cos(t/5.-.5)*.5);
    init = erot(init, vec3(0,0,1), cos(t/5.)*.5);
    cam = erot(cam, vec3(1,0,0), (asin(sin(t/2.))*-1.)*.2);
    init.z -= 2.; 
    
    vec3 p = init;
    bool hit = false;
    float dist;
    float glo = 0.;
    float glo2 = 0.;
    for (int i = 0; i < 100 && !hit; i++) {
        dist = scene(p);
        hit = dist*dist < 1e-6;
        p+=dist*cam;
        glo += 100./(1.+l1*2000.)*dist;
        glo2 += 200./(1.+torch*5000.)*dist;
    }
    float dbgg = dbg;
    bool tc = torch == dist;
    bool pd = dist == pedistal;
    vec3 lid = id;
    vec3 n = norm(p);
    vec3 r = reflect(cam,n);
    
    vec3 lookat = normalize(init-lid);
    vec3 lookat2 = normalize(lpos-lid);
    lookat = normalize(mix(lookat,lookat2,smoothstep(8.,5.,distance(lid,lpos))));

    vec3 eycrd = urot(n, lookat, vec3(1,0,0), vec3(0,0,1));
    float ey = eycrd.x;
    float garble = atan(eycrd.y,eycrd.z) + dot(sin(eycrd.yz*11.),sin(eycrd.yz*25.))*.05;
    garble = (sin(garble*50.)+sin(garble*131.))*.2+.9;
    
    vec3 ldir = normalize(lpos-p);
    float ao = smoothstep(12.,6.,length(p));
    float nd = dot(ldir,n)*.5+.5;
    float rd = max(0.,dot(ldir,n));
    float fres = 1.-abs(dot(cam,n))*.98;
    vec3 ecol = vec3(0.3,0.4,0.7);
    float sd1 = hash(hash(lid.x,lid.y),lid.z);
    float sd2 = hash(sd1,sd1);
    ecol = erot(ecol, vec3(0,1,0), sd2);
    ecol = erot(ecol, vec3(1,0,0), sd1*.4)*garble;
    float atten = 3./pow(distance(p,lpos),1.4) + .4/pow(length(vec2(length(p.xy)-1.8, p.z+3.2)), 1.2);
    
    vec3 dcol = mix(vec3(.9),ecol, smoothstep(.8,.9,ey));
    dcol = mix(dcol,vec3(0), smoothstep(.95,.98,ey));
    
    if (pd) dcol = vec3(0.5,0.4,0.3)*(dbgg*.25+.75);
    vec3 col = (dcol*nd + pow(smoothstep(.7,1.,rd),100.)*fres*1.2)*ao*atten;
    if (distance(lpos,p) < .4) {
        col = vec3(0.2,0.5,0.9);
    }
    if (tc) col = vec3(0.8,0.5,0.2);
    fragColor.xyz = (hit ? col : vec3(.1))+glo*glo + glo*vec3(0.2,0.5,0.9) +glo2*glo2 + glo2*vec3(0.8,0.5,0.2);
    fragColor = smoothstep(vec4(0.1),vec4(1.02),sqrt(fragColor));
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}