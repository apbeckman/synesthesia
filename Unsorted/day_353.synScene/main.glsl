vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** Common Code Begins ********

#define T(u) texture(syn_UserImage,(u)/RENDERSIZE.xy)
#define T1(u) texture(syn_UserImage,(u)/RENDERSIZE.xy)
#define T2(u) texture(syn_UserImage,(u)/RENDERSIZE.xy)

#define R (RENDERSIZE.xy)

#define pal(a,b,c,d,e) ((a) + (b)*sin((c)*(d) + (e)))

			//******** BuffA Code Begins ********


vec4 r42(vec2 u){
    return texture(image30,u/RENDERSIZE.xy);
}

const float velScale = 0.001;


vec4 get(vec2 u){
    vec4 g = T(u);
    vec2 muvn = (iMouse.xy - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    g.zw = mix(g.zw, -normalize(muvn - g.xy ),0.02*float(iMouse.z > 0.)*smoothstep(1.,0.,length(g.xy - muvn)*2. ));
    
    g.xy += sin(g.xy*9. + TIME*3. + g.zw*sin(TIME + g.z))*0.0004;
    
    g.xy += g.zw*velScale;
    
    if (abs(g.x)>0.5*R.x/R.y){
        g.z *= -1.;
        g.x += 2.*g.z*velScale;
    }
    if (abs(g.y)>0.9*R.y/R.x){
        g.w *= -1.;
        g.y += 2.*g.w*velScale;
    }
    
    return g;
}

/*

vec4 otherc = other;                \
        otherc.xy *= 1. - 2.*vec2(             \
            float(abs(me.x) > 0.5 && abs(other.x) > 0.5 && sign(other.x) != sign(me.x) ), \
            float(abs(me.y) > 0.5 && abs(other.y) > 0.5 && sign(other.y) != sign(me.y) )  \
        ); \
        

*/

#define check(dir,times) \
    for(float i = 1.; i < times + 1.; i++){ \
        vec4 other = get(U + dir*i);        \
        if(length(uvn - other.xy) < length(uvn-me.xy)){me = other;}\
    }


vec4 renderPassA() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;

    vec2 uv = (U)/R;
    vec2 uvn = (U - R*0.5)/R.y;


    

    C = T(U);

    vec4 me = get(U);
    
    
    float steps = 4.;
    
    check(vec2(0,1)    ,steps);
    check(vec2(0,-1)   ,steps);
    check(vec2(1,0)    ,steps);
    check(vec2(-1,0)   ,steps);
    check(vec2(1,1)    ,steps);
    check(vec2(-1,-1)  ,steps);
    check(vec2(1,-1)   ,steps);
    check(vec2(-1,1)   ,steps);
    
    C = me;
    

    if(FRAMECOUNT < 2){
        C = r42(U*10.)*2. - 1.;
        //C = vec4(sin(U.x/10.),0.0,1.0,1.0);
    }
	return C; 
 } 


			//******** BuffB Code Begins ********

vec2 uvn;
float stepSz = 1.;

float check(vec4 me,vec2 U,vec2 dir, float times){
    for(float i = 1.; i < times + 1.; i++){
        vec4 other = T(U + dir*i*stepSz); 
        if(other != me){return i;}
    }
    return times;
}

vec4 renderPassB() {
	vec4 col = vec4(0.0);
	vec2 U = _xy;

    vec2 uv = U/R;
 
    uvn = (U - R*0.5)/R.y;

    vec4 me = T(U);

    col = T1(U);
    vec4 C = pal(0.,1.,vec4(5,3,4,1.),1.  + TIME*0.3,4. + me.xyyy*1. + me.z - TIME*0.3 - sin(me.w*21.));
    
    float steps = 14. + sin(TIME + uv.x)*4.;
    
    steps *= 1.4;
    #define shade(expression) C *= (steps - expression)/steps;
    shade(check(me, U, vec2(0,1)    ,steps));
    shade(check(me, U, vec2(0,-1)   ,steps));
    shade(check(me, U, vec2(1,0)    ,steps));
    shade(check(me, U, vec2(-1,0)   ,steps));
    shade(check(me, U, vec2(1,1)    ,steps));
    shade(check(me, U, vec2(-1,-1)  ,steps));
    shade(check(me, U, vec2(1,-1)   ,steps));
    shade(check(me, U, vec2(-1,1)   ,steps));    
    
    
    //C = 1. - C;
    
    //C *= 1.-smoothstep(dFdx(uv.x),0., length(uvn - me.xy) - 0.003);
    //C = mix(C,vec4(1.)*0.01,smoothstep(dFdx(uv.x),0., length(uvn - me.xy) - 0.003));
    
    //C = mix(C,1. - col,smoothstep(dFdx(uv.x),0., length(uvn - me.xy) - 0.003  ));
    
    //C = smoothstep(0.,1.,C);
    
    C = pow(abs(C),vec4(0.4545));
    if (FRAMECOUNT > 2)
        col = mix(C,col,0.4);

	return col; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
	vec4 col = vec4(0.0);
	vec2 U = _xy;

    vec2 uv = U/R;
 
    vec2 uvn = (U - R*0.5)/R.y;

    vec4 me = T(U);

    col = T1(U);
    //C = 1. - C;
    
    //C *= 1.-smoothstep(dFdx(uv.x),0., length(uvn - me.xy) - 0.003);
    //C = mix(C,vec4(1.)*0.01,smoothstep(dFdx(uv.x),0., length(uvn - me.xy) - 0.003));
    
    //C = smoothstep(0.,1.,C);
    
    col = mix(
        col,
        mix(vec4(0),1. - col,smoothstep(dFdx(uv.x),0., length(uvn - me.xy) - 0.003  ))
        ,0.03
    );
    
    //col = pow(abs(col),vec4(0.4545));
    
    
    if (FRAMECOUNT < 2)
        col = col - col;
	return col; 
 } 


vec4 renderMainImage() {
	vec4 C = vec4(0.0);
	vec2 U = _xy;

    vec2 uv = U/R;
    vec2 uvn = (U - R*0.5)/R.y;
    C = T(U);
    
    //C = mix(C,4.-  C, T2(U).x);
    vec4 me = T1(U);
	return C; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderPassC();
	}
	if(PASSINDEX == 3){
		return renderMainImage();
	}
}