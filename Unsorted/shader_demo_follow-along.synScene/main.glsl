

vec4 renderMain(void)

{

    float t = TIME*0.05;

    vec3 color= vec3(0.);

    float circle = 0.;

    vec2 uv = _uvc;





    for(float i = 0.;i<41;i++){

    //uv = uv - (_noise(vec3(TIME*0.1 + uv, TIME*0.1)*2.-1.)*0.25);

    circle = length(uv);

    color += smoothstep(0.9, 0.4, circle)*vec3(_pulse(circle, fract(t + 3.33*i), 0.1));

    color += smoothstep(0.2, 0.5, circle)*vec3(_pulse(1.0-circle, fract(t + 3.33*i), 0.1));



    vec2 circle_position = _uvc - circleposition;



    //float center_dist = length(circle_position - i*0.2);

    //center_dist = smoothstep(slider1+0.05, slider1-0.05, center_dist);

 	//color += vec3(center_dist);

    }

    color -= vec3(sin(color));



	return vec4(color, 1.0);



}

