#include "GL\glew.h"
#include "GLFW\glfw3.h"

#include <vector>
#include "math_utils.h"
#include "shader.h"

#define internal static
#define global_variable static

#define Assert(Expression) if(!(Expression)) { *(int *)0 = 0;}
#define ArrayCount(Array) (sizeof(Array)/sizeof((Array)[0]))

struct body
{
	vec2 P;
	vec2 dP;

	float Orientation;
	float AngularSpeed;

	vec2 Force;
	float Torque;

	float Mass, InvMass;
	float Inertia, InvInertia;
	float CoeffOfRestitution;
	float CoeffOfFriction;

	float Width, Height;
};

struct contact
{
	vec2 GlobalA, GlobalB;
	vec2 N;
	vec2 LocalA, LocalB;
	float PenetrationDepth;

	float AccumulatedNormalImpulse;
	float AccumulatedTangentImpulse;
};
struct manifold
{
	body* BodyA;
	body* BodyB;

	uint32_t ContactCount;
	contact Contacts[2];

	manifold(body* A, body* B) { BodyA = A; BodyB = B; ContactCount = 0; }
	void ChooseBestPoints(contact& NewContact);
};

struct world
{
	void AddBody(body &Body);

	std::vector<body> Bodies;
	std::vector<manifold> Manifolds;
};

void world::AddBody(body &Body)
{
	Bodies.push_back(Body);
}

void manifold::ChooseBestPoints(contact &NewContact)
{
	contact TempArray[3] = { Contacts[0], Contacts[1], NewContact };

	float Dist01 = LengthSq(TempArray[0].GlobalA - TempArray[1].GlobalA);
	float Dist02 = LengthSq(TempArray[0].GlobalA - TempArray[2].GlobalA);
	float Dist12 = LengthSq(TempArray[1].GlobalA - TempArray[2].GlobalA);

#if 1
	// NOTE(georgy): This path chooses 2 most distant from each other points
	if(Dist01 > Dist02)
	{
		if(Dist01 > Dist12)
		{
			Contacts[0] = TempArray[0];
			Contacts[1] = TempArray[1];
		}
		else
		{
			Contacts[0] = TempArray[1];
			Contacts[1] = TempArray[2];
		}
	}
	else
	{
		if(Dist02 > Dist12)
		{
			Contacts[0] = TempArray[0];
			Contacts[1] = TempArray[2];
		}
		else
		{
			Contacts[0] = TempArray[1];
			Contacts[1] = TempArray[2];
		}
	}
#else
	// NOTE(georgy): This path chooses the most penetrated point as the first point 
	// 				 and the second point is the most distant from the first
	float Penetration = -FLT_MAX;
	for(uint32_t ContactIndex = 0;
		ContactIndex < ArrayCount(TempArray);
		ContactIndex++)
	{
		contact *Contact = TempArray + ContactIndex;

		if(Contact->PenetrationDepth > Penetration)
		{
			Penetration = Contact->PenetrationDepth;
			Contacts[0] = *Contact;
		}
	}

	float DistanceSq = -FLT_MAX;
	for(uint32_t ContactIndex = 0;
		ContactIndex < ArrayCount(TempArray);
		ContactIndex++)
	{
		contact *Contact = TempArray + ContactIndex;
		float Dist = LengthSq(Contact->GlobalA - Contacts[0].GlobalA);

		if(Dist > DistanceSq)
		{
			DistanceSq = Dist;
			Contacts[1] = *Contact;
		}
	}
#endif
}

internal int32_t
FindManifold(world* World, body* A, body* B)
{
	int32_t Result = -1;

	for (uint32_t i = 0; i < World->Manifolds.size(); i++)
	{
		manifold* Manifold = &World->Manifolds[i];
		if ((A == Manifold->BodyA) && (B == Manifold->BodyB))
		{
			Result = i;
			break;
		}
	}

	return(Result);
}


internal manifold *
GetOrCreateManifold(world *World, body *BodyA, body *BodyB)
{
	manifold *Result = 0;

	int32_t ManifoldIndex = FindManifold(World, BodyA, BodyB);
	if(ManifoldIndex != -1)
	{
		Result = &World->Manifolds[ManifoldIndex];
	}
	else
	{
		manifold NewManifold(BodyA, BodyB);
		World->Manifolds.push_back(NewManifold);

		Result = &World->Manifolds[World->Manifolds.size() - 1];
	}

	return(Result);
}

global_variable world World;

internal void
RenderBody(body& Body, shader Shader)
{
	mat4 Model = Translation(vec3(Body.P, 0.0f)) *
				 Rotation(Body.Orientation, vec3(0.0f, 0.0f, 1.0f)) *
				 Scaling(vec3(Body.Width, Body.Height, 1.0));
	Shader.SetMat4("Model", Model);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}

static uint32_t WindowWidth;
static uint32_t WindowHeight;
static void
GLFWFramebufferSizeCallback(GLFWwindow* Window, int NewWidth, int NewHeight)
{
	WindowWidth = NewWidth;
	WindowHeight = NewHeight;
	glViewport(0, 0, NewWidth, NewHeight);
}

internal void
GLFWKeyCallback(GLFWwindow* Window, int Key, int Scancode, int Action, int Mods)
{}

internal vec2
Support(body& Body, vec2 D)
{
	float Rad = Radians(Body.Orientation);
	vec2 XAxis = vec2(cosf(Rad), sinf(Rad));
	vec2 YAxis = Perp(XAxis);
	float HalfWidth = 0.5f * Body.Width;
	float HalfHeight = 0.5f * Body.Height;

	vec2 Points[4] =
	{
		Body.P - HalfWidth * XAxis + HalfHeight * YAxis,
		Body.P - HalfWidth * XAxis - HalfHeight * YAxis,
		Body.P + HalfWidth * XAxis - HalfHeight * YAxis,
		Body.P + HalfWidth * XAxis + HalfHeight * YAxis,
	};

	vec2 Result = Points[0];
	float Max = Dot(D, Points[0]);
	for (uint32_t i = 1; i < 4; i++)
	{
		float Value = Dot(D, Points[i]);
		if (Value > Max)
		{
			Max = Value;
			Result = Points[i];
		}
	}

	return(Result);
}

internal bool
DoSimplex(vec2* Simplex, vec2* SimplexA, vec2* SimplexB, uint32_t& SimplexCount, vec2& D)
{
	switch (SimplexCount)
	{
		case 2:
		{
			vec2 A = Simplex[1];
			vec2 B = Simplex[0];
			vec2 AB = B - A;
			vec2 AO = vec2(0.0f, 0.0f) - A;

			if(Dot(AB, AO) > 0.0)
			{
				if(Cross2D(AB, AO) > 0.0f)
				{
					D = Perp(AB);
				}
				else
				{
					D = -Perp(AB);
				}
			}
			else
			{
				SimplexCount = 1;
				Simplex[0] = Simplex[1];
				SimplexA[0] = SimplexA[1];
				SimplexB[0] = SimplexB[1];
			}
		} break;

		case 3:
		{
			vec2 A = Simplex[2];
			vec2 B = Simplex[1];
			vec2 C = Simplex[0];
			vec2 AB = B - A;
			vec2 AC = C - A;
			vec2 AO = vec2(0.0f, 0.0f) - A;

			// NOTE(georgy): I want triangle to be CCW
			if (Cross2D(AB, AC) < 0.0f)
			{
				vec2 Temp = Simplex[0];
				Simplex[0] = Simplex[1];
				Simplex[1] = Temp;

				Temp = SimplexA[0];
				SimplexA[0] = SimplexA[1];
				SimplexA[1] = Temp;

				Temp = SimplexB[0];
				SimplexB[0] = SimplexB[1];
				SimplexB[1] = Temp;

				B = Simplex[1];
				C = Simplex[0];
				AB = B - A;
				AC = C - A;
			}
			Assert(Cross2D(AB, AC) > 0.0f);

			vec2 ABNormal = -Perp(AB);
			vec2 ACNormal = Perp(AC);

			if (Dot(ACNormal, AO) > 0.0f)
			{
				if (Dot(AC, AO) > 0.0f)
				{
					SimplexCount = 2;
					Simplex[1] = A;
					SimplexA[1] = SimplexA[2];
					SimplexB[1] = SimplexB[2];
					D = ACNormal;
				}
				else if (Dot(AB, AO) > 0.0f)
				{
					SimplexCount = 2;
					Simplex[0] = B;
					Simplex[1] = A;
					SimplexA[0] = SimplexA[1];
					SimplexA[1] = SimplexA[2];
					SimplexB[0] = SimplexB[1];
					SimplexB[1] = SimplexB[2];
					D = ABNormal;
				}
				else
				{
					SimplexCount = 1;
					Simplex[0] = A;
					SimplexA[0] = SimplexA[2];
					SimplexB[0] = SimplexA[2];
				}
			}
			else
			{
				if (Dot(ABNormal, AO) > 0.0f)
				{
					if (Dot(AB, AO) > 0.0f)
					{
						SimplexCount = 2;
						Simplex[0] = B;
						Simplex[1] = A;
						SimplexA[0] = SimplexA[1];
						SimplexA[1] = SimplexA[2];
						SimplexB[0] = SimplexB[1];
						SimplexB[1] = SimplexB[2];
						D = ABNormal;
					}
					else
					{
						SimplexCount = 1;
						Simplex[0] = A;
						SimplexA[0] = SimplexA[2];
						SimplexB[0] = SimplexA[2];
					}
				}
				else
				{
					return(true);
				}
			}
		} break;
	}

	return(false);
}

struct gjk_result
{
	bool Intersection;

	union
	{
		// NOTE(georgy): If no intersection, returns closest points
		struct
		{
			vec2 ClosestPA, ClosestPB;
		};
		// NOTE(georgy): If intersection, returns the last simplex
		struct 
		{
			vec2 Simplex[3];
			vec2 SimplexA[3];
			vec2 SimplexB[3];
		};
	};
};
internal gjk_result
GJK(body& BodyA, body& BodyB)
{
	gjk_result Result = {};

	uint32_t SimplexCount = 1;
	Result.SimplexB[0] = Support(BodyB, vec2(1.0f, 0.0f));
	Result.SimplexA[0] = Support(BodyA, -vec2(1.0f, 0.0f));
	Result.Simplex[0] = Result.SimplexB[0] - Result.SimplexA[0];
	vec2 D = -Result.Simplex[0];
	while (true)
	{
		vec2 Sb = Support(BodyB, D);
		vec2 Sa = Support(BodyA, -D);
		vec2 S = Sb - Sa;

		// NOTE(georgy): Checks if we already have the farthest point in this direction
		for (uint32_t i = 0; i < SimplexCount; i++)
		{
			if (VectorsAreEqual(S, Result.Simplex[i]))
			{
				if (SimplexCount == 1)
				{
					Result.ClosestPA = Result.SimplexA[0];
					Result.ClosestPB = Result.SimplexB[0];
				}
				else if (SimplexCount == 2)
				{
					vec2 A = Result.Simplex[1];
					vec2 B = Result.Simplex[0];
					vec2 AB = B - A;
					vec2 AO = vec2(0.0f, 0.0f) - A;

					float t = Dot(AO, AB) / Dot(AB, AB);
					Result.ClosestPA = Result.SimplexA[1] + t * (Result.SimplexA[0] - Result.SimplexA[1]);
					Result.ClosestPB = Result.SimplexB[1] + t * (Result.SimplexB[0] - Result.SimplexB[1]);
				}

				return(Result);
			}
		}

		// NOTE(georgy): Checks for degenerate triangle
		if (SimplexCount == 2)
		{
			vec2 AB = Result.Simplex[1] - S;
			vec2 AC = Result.Simplex[0] - S;
			if (fabs(Cross2D(AB, AC)) <= 5.0f * Epsilon)
			{
				vec2 A = Result.Simplex[1];
				vec2 B = Result.Simplex[0];
				vec2 AB = B - A;
				vec2 AO = vec2(0.0f, 0.0f) - A;

				float t = Dot(AO, AB) / Dot(AB, AB);
				Result.ClosestPA = Result.SimplexA[1] + t * (Result.SimplexA[0] - Result.SimplexA[1]);
				Result.ClosestPB = Result.SimplexB[1] + t * (Result.SimplexB[0] - Result.SimplexB[1]);

				return(Result);
			}
		}
		Result.SimplexA[SimplexCount] = Sa;
		Result.SimplexB[SimplexCount] = Sb;
		Result.Simplex[SimplexCount++] = S;
		if (DoSimplex(Result.Simplex, Result.SimplexA, Result.SimplexB, SimplexCount, D))
		{
			Result.Intersection = true;
			return(Result);
		}
	}
}

internal void
EPA(gjk_result* GJKInfo, manifold* Manifold)
{
	const float Tolerance = 0.01f;
	const vec2 O = vec2(0.0f, 0.0);

	std::vector<vec2> Polygon;
	Polygon.push_back(GJKInfo->Simplex[2]);
	Polygon.push_back(GJKInfo->Simplex[1]);
	Polygon.push_back(GJKInfo->Simplex[0]);
	std::vector<vec2> PolygonA;
	PolygonA.push_back(GJKInfo->SimplexA[2]);
	PolygonA.push_back(GJKInfo->SimplexA[1]);
	PolygonA.push_back(GJKInfo->SimplexA[0]);
	std::vector<vec2> PolygonB;
	PolygonB.push_back(GJKInfo->SimplexB[2]);
	PolygonB.push_back(GJKInfo->SimplexB[1]);
	PolygonB.push_back(GJKInfo->SimplexB[0]);
	while (true)
	{
		float ClosestDistance = FLT_MAX;
		vec2 ClosestN = vec2(0.0f, 0.0f);
		int32_t ClosestEdgeIndex = -1;
		for (uint32_t I0 = 0, I1 = Polygon.size() - 1;
			I0 < Polygon.size();
			I1 = I0, I0++)
		{
			vec2 Edge = Polygon[I0] - Polygon[I1];
			vec2 N = Normalize(-Perp(Edge));
			float Distance = Dot(N, Polygon[I0] - O);
			if(Distance < ClosestDistance)
			{
				ClosestDistance = Distance;
				ClosestN = N;
				ClosestEdgeIndex = I1;
			}
		}

		vec2 bP = Support(*Manifold->BodyB, ClosestN);
		vec2 aP = Support(*Manifold->BodyA, -ClosestN);
		vec2 P = bP - aP;
		float Distance = Dot(P, ClosestN);
		if(Absolute(Distance - ClosestDistance) < Tolerance)
		{
			vec2 ProjectedO = O + ClosestDistance*ClosestN;

			uint32_t NextVertexIndex = (ClosestEdgeIndex + 1) % Polygon.size();
			vec2 Edge = Polygon[NextVertexIndex] - Polygon[ClosestEdgeIndex];
			float t = Dot(Edge, ProjectedO - Polygon[ClosestEdgeIndex]) / LengthSq(Edge);

			contact NewContact;
			NewContact.GlobalA = PolygonA[ClosestEdgeIndex] + t*(PolygonA[NextVertexIndex] - PolygonA[ClosestEdgeIndex]);
			NewContact.GlobalB = PolygonB[ClosestEdgeIndex] + t*(PolygonB[NextVertexIndex] - PolygonB[ClosestEdgeIndex]);
			mat2 InvOrientationA = Rotation2x2(-Manifold->BodyA->Orientation);
			NewContact.LocalA = InvOrientationA*(NewContact.GlobalA - Manifold->BodyA->P);
			mat2 InvOrientationB = Rotation2x2(-Manifold->BodyB->Orientation);
			NewContact.LocalB = InvOrientationB*(NewContact.GlobalB - Manifold->BodyB->P);
			NewContact.N = ClosestN;
			NewContact.PenetrationDepth = ClosestDistance;
			NewContact.AccumulatedNormalImpulse = 0.0f;
			NewContact.AccumulatedTangentImpulse = 0.0f;

			bool FarEnough = true;
			for(uint32_t ContactIndex = 0;
				ContactIndex < Manifold->ContactCount;
				ContactIndex++)
			{
				contact *Contact = Manifold->Contacts + ContactIndex;
				vec2 rA = NewContact.GlobalA - Contact->GlobalA;
				vec2 rB = NewContact.GlobalB - Contact->GlobalB;

				const float PersistentThresholdSq = Square(0.1f);
				bool rATooClose = (LengthSq(rA) <= PersistentThresholdSq);
				bool rBTooClose = (LengthSq(rB) <= PersistentThresholdSq);

				if(rATooClose || rBTooClose)
				{
					Contact->GlobalA = NewContact.GlobalA;
					Contact->GlobalB = NewContact.GlobalB;
					Contact->LocalA = NewContact.LocalA;
					Contact->LocalB = NewContact.LocalB;
					Contact->N = NewContact.N;
					Contact->PenetrationDepth = NewContact.PenetrationDepth;

					FarEnough = false;
					break;
				}
			}
			if(FarEnough)
			{
				if(Manifold->ContactCount == 2)
				{
					Manifold->ChooseBestPoints(NewContact);
				}
				else
				{
					Manifold->Contacts[Manifold->ContactCount++] = NewContact;
				}
			}

			break;
		}
		else
		{
			std::vector<vec2>::iterator Iter = Polygon.begin() + (ClosestEdgeIndex + 1);
			Polygon.insert(Iter, P);

			Iter = PolygonA.begin() + (ClosestEdgeIndex + 1);
			PolygonA.insert(Iter, aP);

			Iter = PolygonB.begin() + (ClosestEdgeIndex + 1);
			PolygonB.insert(Iter, bP);
		}
	}
}

int main()
{
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SAMPLES, 4);

	WindowWidth = 900; WindowHeight = 540;
	GLFWwindow* Window = glfwCreateWindow(WindowWidth, WindowHeight, "Physics2D", 0, 0);
	glfwMakeContextCurrent(Window);
	glfwSwapInterval(1);
	glfwSetFramebufferSizeCallback(Window, GLFWFramebufferSizeCallback);
	glfwSetKeyCallback(Window, GLFWKeyCallback);

	GLFWmonitor* Monitor = glfwGetPrimaryMonitor();
	const GLFWvidmode* VidMode = glfwGetVideoMode(Monitor);
	float TargetSecondsPerFrame = 1.0f / VidMode->refreshRate;

	glewInit();

	glViewport(0, 0, WindowWidth, WindowHeight);
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	glClearColor(0.2f, 0.4f, 0.8f, 1.0f);

	float QuadVertices[] =
	{
		-0.5f, 0.5f,
		-0.5f, -0.5f,
		0.5f, 0.5f,
		0.5f, -0.5f
	};
	GLuint QuadVAO, QuadVBO;
	glGenVertexArrays(1, &QuadVAO);
	glGenBuffers(1, &QuadVBO);
	glBindVertexArray(QuadVAO);
	glBindBuffer(GL_ARRAY_BUFFER, QuadVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(QuadVertices), QuadVertices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
	glBindVertexArray(0);

#define FIRST_DEMO 0
#define SECOND_DEMO 1
#define THIRD_DEMO 0
#if FIRST_DEMO
	for(int32_t Y = -8; Y <= 5; Y++)
	{
		body Body = {};
		Body.P = vec2(0.0f, (float)Y);
		Body.Width = 1.0f;
		Body.Height = 1.0f;
		Body.Mass = 10.0f;
		Body.InvMass = 1.0f / Body.Mass;
		Body.Inertia = (1.0f / 12.0f) * Body.Mass * (Body.Width * Body.Width + Body.Height * Body.Height);
		Body.InvInertia = 1.0f / Body.Inertia;
		Body.Orientation = 0.0f;
		Body.CoeffOfRestitution = 0.0f;
		Body.CoeffOfFriction = 0.2f;
		World.AddBody(Body);
	}
#elif SECOND_DEMO
	for(int32_t Y = -2; Y <= 6; Y += 2)
	{
		for(int32_t X = -6; X <= 6; X += 2)
		{
			body Body = {};
			Body.P = vec2((float)X, (float)Y);
			Body.Width = 1.0f;
			Body.Height = 1.0f;
			Body.Orientation = 18.8f;
			Body.Mass = 10.0f;
			Body.InvMass = 1.0f / Body.Mass;
			Body.Inertia = (1.0f / 12.0f) * Body.Mass * (Body.Width * Body.Width + Body.Height * Body.Height);
			Body.InvInertia = 1.0f / Body.Inertia;
			Body.CoeffOfRestitution = 0.2f;
			Body.CoeffOfFriction = 0.2f;
			World.AddBody(Body);
		}
	}
#elif THIRD_DEMO
	vec2 Left = vec2(-5.0f, -7.0f);
	for(int32_t i = 0; i < 10; i++)
	{
		vec2 P = Left;

		for(int32_t j = i; j < 10; j++)
		{
			body Body = {};
			Body.P = P;
			Body.Width = 1.0f;
			Body.Height = 1.0f;
			Body.Mass = 10.0f;
			Body.InvMass = 1.0f / Body.Mass;
			Body.Inertia = (1.0f/12.0f)*Body.Mass*(Body.Width*Body.Width + Body.Height*Body.Height);
			Body.InvInertia = 1.0f / Body.Inertia;
			Body.CoeffOfRestitution = 0.0f;
			Body.CoeffOfFriction = 1.0f;
			World.AddBody(Body);

			P += vec2(1.0f, 0.0f);
		}

		Left += vec2(0.5, 2.0f);
	}
#endif

	body Body = {};
	Body.P = vec2(0.0f, -12.5f);
	Body.Mass = FLT_MAX;
	Body.InvMass = 0.0f;
	Body.Inertia = FLT_MAX;
	Body.InvInertia = 0.0f;
	Body.Width = 100.0f;
	Body.Height = 8.0f;
	Body.Orientation = 0.0f;
	Body.CoeffOfRestitution = 0.0f;
	Body.CoeffOfFriction = 0.2f;
	World.AddBody(Body);

	shader DefaultShader("shaders\\vertex.vert", "shaders\\fragment.frag");
	float dt = TargetSecondsPerFrame;
	float LastTime = (float)glfwGetTime();
	while(!glfwWindowShouldClose(Window))
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		mat4 Projection;
		float AspectRatio = (float)WindowWidth / (float)WindowHeight;
		if (AspectRatio >= 1.0f)
		{
			Projection = Ortho(-10.0f, 10.0f, -10.0f * AspectRatio, 10.0f * AspectRatio, 0.0f, 50.0f);
		}
		else
		{
			Projection = Ortho(-10.0f / AspectRatio, 10.0f / AspectRatio, -10.0f, 10.0f, 0.0f, 50.0f);
		}

		for(uint32_t i = 0; i < World.Bodies.size(); i++)
		{
			body *Body = &World.Bodies[i];

			if(Body->InvMass != 0.0f)
			{
				Body->dP += dt * (vec2(0.0f, -9.8f) + Body->Force * Body->InvMass);
				Body->AngularSpeed += dt * (Body->Torque * Body->InvInertia);
			}
		}	

		for(uint32_t i = 0; i < World.Manifolds.size(); i++)
		{
			manifold *Manifold = &World.Manifolds[i];

			for(uint32_t j = 0; j < Manifold->ContactCount;)
			{
				contact *Contact = Manifold->Contacts + j;
				body *BodyA = Manifold->BodyA;
				body *BodyB = Manifold->BodyB;
				const float PersistentThreshold = 0.1f;
				const float PersistentThresholdSq = Square(PersistentThreshold);

				vec2 NewGlobalA = BodyA->P + Rotation2x2(BodyA->Orientation)*Contact->LocalA;
				vec2 NewGlobalB = BodyB->P + Rotation2x2(BodyB->Orientation)*Contact->LocalB;
				Contact->GlobalA = NewGlobalA;
				Contact->GlobalB = NewGlobalB;

				vec2 rAB = (NewGlobalB - NewGlobalA);
				if(Dot(rAB, Contact->N) <= -PersistentThreshold)
				{
					Manifold->Contacts[j] = Manifold->Contacts[--Manifold->ContactCount];
				}
				else
				{
					vec2 ProjectedP = NewGlobalA + Contact->N*Dot(rAB, Contact->N);
					float DistSq = LengthSq(NewGlobalB - ProjectedP);
					if(DistSq > PersistentThresholdSq)
					{
						Manifold->Contacts[j] = Manifold->Contacts[--Manifold->ContactCount];
					}
					else
					{
						j++;
					}
				}
			}
		}
		
		for (uint32_t i = 0; i < World.Bodies.size() - 1; i++)
		{
			body* BodyA = &World.Bodies[i];
			for (uint32_t j = i + 1; j < World.Bodies.size(); j++)
			{
				body* BodyB = &World.Bodies[j];

				gjk_result GJKResult = GJK(*BodyA, *BodyB);
				if(GJKResult.Intersection)
				{
					manifold *Manifold = GetOrCreateManifold(&World, BodyA, BodyB);
					EPA(&GJKResult, Manifold);
				}
				else
				{
					float Margin = 0.01f;
					float LengthBetweenClosest = Length(GJKResult.ClosestPA - GJKResult.ClosestPB);
					if(LengthBetweenClosest < 2.0f*Margin)
					{
						manifold *Manifold = GetOrCreateManifold(&World, BodyA, BodyB);

						contact NewContact;
						NewContact.GlobalA = GJKResult.ClosestPA;
						NewContact.GlobalB = GJKResult.ClosestPB;
						mat2 InvOrientationA = Rotation2x2(-Manifold->BodyA->Orientation);
						NewContact.LocalA = InvOrientationA*(NewContact.GlobalA - Manifold->BodyA->P);
						mat2 InvOrientationB = Rotation2x2(-Manifold->BodyB->Orientation);
						NewContact.LocalB = InvOrientationB*(NewContact.GlobalB - Manifold->BodyB->P);
						NewContact.N = NOZ(GJKResult.ClosestPA - GJKResult.ClosestPB);
						NewContact.PenetrationDepth = Max(Margin - LengthBetweenClosest, 0.0f);
						NewContact.AccumulatedNormalImpulse = 0.0f;
						NewContact.AccumulatedTangentImpulse = 0.0f;

						bool FarEnough = true;
						for(uint32_t ContactIndex = 0;
							ContactIndex < Manifold->ContactCount;
							ContactIndex++)
						{
							contact *Contact = Manifold->Contacts + ContactIndex;
							vec2 rA = NewContact.GlobalA - Contact->GlobalA;
							vec2 rB = NewContact.GlobalB - Contact->GlobalB;

							const float PersistentThresholdSq = Square(0.1f);
							bool rATooClose = (LengthSq(rA) <= PersistentThresholdSq);
							bool rBTooClose = (LengthSq(rB) <= PersistentThresholdSq);

							if(rATooClose || rBTooClose)
							{
								Contact->GlobalA = NewContact.GlobalA;
								Contact->GlobalB = NewContact.GlobalB;
								Contact->LocalA = NewContact.LocalA;
								Contact->LocalB = NewContact.LocalB;
								Contact->N = NewContact.N;
								Contact->PenetrationDepth = NewContact.PenetrationDepth;

								FarEnough = false;
								break;
							}
						}
						if(FarEnough)
						{
							if(Manifold->ContactCount == 2)
							{
								Manifold->ChooseBestPoints(NewContact);
							}
							else
							{
								Manifold->Contacts[Manifold->ContactCount++] = NewContact;
							}
						}
					}
					else
					{
						for(uint32_t i = 0; i < World.Manifolds.size(); i++)
						{
							manifold *Manifold = &World.Manifolds[i];
							if((BodyA == Manifold->BodyA) && (BodyB == Manifold->BodyB))
							{
								World.Manifolds.erase(World.Manifolds.begin() + i);
							}
						}
					}
				}
			}
		}

		for(uint32_t ManifoldIndex = 0; 
			ManifoldIndex < World.Manifolds.size(); 
			ManifoldIndex++)
		{
			manifold* Manifold = &World.Manifolds[ManifoldIndex];
			body* BodyA = Manifold->BodyA;
			body* BodyB = Manifold->BodyB;

			for(uint32_t ContactIndex = 0;
				ContactIndex < Manifold->ContactCount;
				ContactIndex++)
			{
				contact *Contact = Manifold->Contacts + ContactIndex;
				vec2 N = Contact->N;
				vec2 Tangent = Perp(N);
				vec2 Impulse = Contact->AccumulatedNormalImpulse*N + Contact->AccumulatedTangentImpulse*Tangent;

				vec2 aP = Contact->GlobalA - BodyA->P;
				vec2 bP = Contact->GlobalB - BodyB->P;
				BodyA->dP += Impulse * BodyA->InvMass;
				BodyA->AngularSpeed += BodyA->InvInertia * Dot(Perp(aP), Impulse);
				BodyB->dP += -Impulse * BodyB->InvMass;
				BodyB->AngularSpeed += BodyB->InvInertia * Dot(Perp(bP), -Impulse);
			}
		}
		
		for(uint32_t Iter = 0;
			Iter < 20;
			Iter++)
		{
			for(uint32_t ManifoldIndex = 0; 
				ManifoldIndex < World.Manifolds.size(); 
				ManifoldIndex++)
			{
				manifold* Manifold = &World.Manifolds[ManifoldIndex];
				body* BodyA = Manifold->BodyA;
				body* BodyB = Manifold->BodyB;

				float CoeffOfRestitution = Min(BodyA->CoeffOfRestitution, BodyB->CoeffOfRestitution);
				float CoeffOfFriction = sqrt(BodyA->CoeffOfFriction * BodyB->CoeffOfFriction);

				for(uint32_t ContactIndex = 0;
					ContactIndex < Manifold->ContactCount;
					ContactIndex++)
				{
					contact *Contact = Manifold->Contacts + ContactIndex;
					vec2 N = Contact->N;

					float Slop = 0.01f;
					float BiasFactor = 0.2f;
					float BiasVelocity = (BiasFactor / dt) * Max(0.0f, Contact->PenetrationDepth - Slop);

					vec2 aP = Contact->GlobalA - BodyA->P;
					vec2 VelocityA = BodyA->dP + BodyA->AngularSpeed * Perp(aP);

					vec2 bP = Contact->GlobalB - BodyB->P;
					vec2 VelocityB = BodyB->dP + BodyB->AngularSpeed * Perp(bP);

					vec2 RelativeVel = VelocityA - VelocityB;
					float ImpulseNom = -(1.0f + CoeffOfRestitution) * Dot(RelativeVel, N) + BiasVelocity;
					float ImpulseDenom = (BodyA->InvMass + BodyB->InvMass) + 
										  Square(Dot(Perp(aP), N))*BodyA->InvInertia + Square(Dot(Perp(bP), N))*BodyB->InvInertia;
					float ImpulseMagnitude = ImpulseNom / ImpulseDenom;
					
					float TempImpulse = Contact->AccumulatedNormalImpulse;
					Contact->AccumulatedNormalImpulse = Max(Contact->AccumulatedNormalImpulse + ImpulseMagnitude, 0.0f);
					ImpulseMagnitude = Contact->AccumulatedNormalImpulse - TempImpulse;

					BodyB->dP += -ImpulseMagnitude * BodyB->InvMass * N;
					BodyB->AngularSpeed += BodyB->InvInertia * Dot(Perp(bP), -ImpulseMagnitude * N);

					BodyA->dP += ImpulseMagnitude * BodyA->InvMass * N;
					BodyA->AngularSpeed += BodyA->InvInertia * Dot(Perp(aP), ImpulseMagnitude * N);

					vec2 Tangent = Perp(N);
					RelativeVel = (BodyA->dP + BodyA->AngularSpeed*Perp(aP)) - (BodyB->dP + BodyB->AngularSpeed*Perp(bP));
					float ImpulseNomT = -(1.0f + CoeffOfRestitution) * Dot(RelativeVel, Tangent);
					float ImpulseDenomT = (BodyA->InvMass + BodyB->InvMass) + 
					 					   Square(Dot(Perp(aP), Tangent))*BodyA->InvInertia + Square(Dot(Perp(bP), Tangent))*BodyB->InvInertia;
					float ImpulseMagnitudeT = ImpulseNomT / ImpulseDenomT;

					float TempImpulseT = Contact->AccumulatedTangentImpulse;
					Contact->AccumulatedTangentImpulse = Clamp(Contact->AccumulatedTangentImpulse + ImpulseMagnitudeT,
															   -CoeffOfFriction*ImpulseMagnitude, CoeffOfFriction*ImpulseMagnitude);
					ImpulseMagnitudeT = Contact->AccumulatedTangentImpulse - TempImpulseT;

					BodyB->dP += -ImpulseMagnitudeT * BodyB->InvMass * Tangent;
					BodyB->AngularSpeed += BodyB->InvInertia * Dot(Perp(bP), -ImpulseMagnitudeT * Tangent);

					BodyA->dP += ImpulseMagnitudeT * BodyA->InvMass * Tangent;
					BodyA->AngularSpeed += BodyA->InvInertia * Dot(Perp(aP), ImpulseMagnitudeT * Tangent);
				}
			}
		}

		for (uint32_t i = 0; i < World.Bodies.size(); i++)
		{
			body* Body = &World.Bodies[i];

			if (Body->InvMass != 0.0f)
			{
				Body->P += dt * Body->dP;
				Body->Orientation += dt * Degrees(Body->AngularSpeed);
			}

			Body->Force = vec2(0.0f, 0.0f);
			Body->Torque = 0.0f;
		}

		DefaultShader.Use();
		DefaultShader.SetMat4("Projection", Projection);
		glBindVertexArray(QuadVAO);
		for(uint32_t ManifoldIndex = 0; 
			ManifoldIndex < World.Manifolds.size(); 
			ManifoldIndex++)
		{
			manifold* Manifold = &World.Manifolds[ManifoldIndex];
			body* BodyA = Manifold->BodyA;
			body* BodyB = Manifold->BodyB;

			for(uint32_t ContactIndex = 0;
				ContactIndex < Manifold->ContactCount;
				ContactIndex++)
			{
				contact *Contact = Manifold->Contacts + ContactIndex;

				mat4 Model = Translation(vec3(Contact->GlobalA, 0.0f)) *
							 Scaling(vec3(0.3f, 0.3f, 1.0));
				DefaultShader.SetMat4("Model", Model);
				DefaultShader.SetVec3("Color", vec3(0.0f, 1.0f, 0.0f));
				glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

				Model = Translation(vec3(Contact->GlobalB, 0.0f)) *
						Scaling(vec3(0.3f, 0.3f, 1.0));
				DefaultShader.SetMat4("Model", Model);
				DefaultShader.SetVec3("Color", vec3(1.0f, 0.0f, 0.0f));
				glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
			}
		}

		DefaultShader.SetVec3("Color", vec3(1.0f, 1.0f, 1.0f));
		for (int i = 0; i < World.Bodies.size(); i++)
		{
			RenderBody(World.Bodies[i], DefaultShader);
		}
		glBindVertexArray(0);

		glfwPollEvents();
		glfwSwapBuffers(Window);

		if (((float)glfwGetTime() - LastTime) < TargetSecondsPerFrame)
		{
			while (((float)glfwGetTime() - LastTime) < TargetSecondsPerFrame)
			{}
		}
		std::cout << glfwGetTime() - LastTime << std::endl;
		LastTime = glfwGetTime();
	}

	return(0);
}