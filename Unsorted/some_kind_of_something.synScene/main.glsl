

			//******** Common Code Begins ********


#define R RENDERSIZE.xy

#define T(u) texture(iChannel0,(u)/R)
#define T1(u) texture(iChannel1,(u)/R)
#define T2(u) texture(iChannel2,(u)/R)
#define T3(u) texture(iChannel3,(u)/R)


#define Neighbors vec4 n = T(U + vec2(0,1)); vec4 s = T(U - vec2(0,1)); vec4 e = T(U + vec2(1,0)); vec4 w = T(U - vec2(1,0)); vec4 m = 0.25*(n + w + e + s);




			//******** BuffA Code Begins ********

mat2 rot(float a){
    float s = sin(a), c = cos(a);
    return mat2(c,s,-s,c);
}


vec4 renderPassA() {
	vec4 C = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy;
    
    if(!bool(FRAMECOUNT==20))
    {
        C = vec4(1);
        if(abs(uv.x-.5) > 0.45 || abs(uv.y-.5) > .45)
        {
            C = vec4(0.);
        }
        uv = abs(uv-.5);
		//uv -= .5;
        float a = .1*sin(smoothTime*0.2);
        float b = .1*sin(smoothTime*.546);

        uv *= rot(a);
        uv /= .9+b;

        uv += .5;

        //C *= .5+.5*sin(vec4(1,2,3,4)*TIME*84.+vec4(34,23,562,0));

        C *= (1.-texture(BuffA,uv))*1.5;
        
    }
    else
        C = texture(BuffA,uv);
    
	return C; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    fragCoord += 2.*sin(fragCoord.yx*.01+TIME);
    
    fragCoord += 3.*sin(TIME*94.*vec2(3434,253));
    if(!bool((FRAMECOUNT+int(fragCoord.x))==10) )
    	for(int i = 0; i < 4; i++)
 	   		fragColor[i] = texture(BuffA,(fragCoord+2.*float(i))/RENDERSIZE.xy)[i];
    else
        discard;
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}