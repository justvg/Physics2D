#version 330 core
layout (location = 0) in vec2 aP;
 
uniform mat4 Projection = mat4(1.0);
uniform mat4 View = mat4(1.0);
uniform mat4 Model = mat4(1.0);

void main()
{
    gl_Position = Projection * View * Model * vec4(aP, 0.0, 1.0);
}