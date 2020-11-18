// OBJ_Loader.h - A Single Header OBJ Model Loader

#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <sstream>
#include <locale>

#include <math.h>

// Print progress to console while loading (large models)
// #define OBJL_CONSOLE_OUTPUT

// Namespace: OBJL
//
// Description: The namespace that holds eveyrthing that
//	is needed and used for the OBJ Model Loader
namespace objl
{

namespace
{
auto readFloat = [](const std::string& s)
{
  float v;
  std::istringstream istr(s);
  istr.imbue(std::locale("C"));
  istr >> v;
  return v;
};
}

// Structure: Vector2
//
// Description: A 2D Vector that Holds Positional Data
struct Vector2
{
  // Default Constructor
  Vector2()
  {
    X = 0.0f;
    Y = 0.0f;
  }
  // Variable Set Constructor
  Vector2(float X_, float Y_)
  {
    X = X_;
    Y = Y_;
  }
  // Bool Equals Operator Overload
  bool operator==(const Vector2& other) const
  {
    return (this->X == other.X && this->Y == other.Y);
  }
  // Bool Not Equals Operator Overload
  bool operator!=(const Vector2& other) const
  {
    return !(this->X == other.X && this->Y == other.Y);
  }
  // Addition Operator Overload
  Vector2 operator+(const Vector2& right) const
  {
    return Vector2(this->X + right.X, this->Y + right.Y);
  }
  // Subtraction Operator Overload
  Vector2 operator-(const Vector2& right) const
  {
    return Vector2(this->X - right.X, this->Y - right.Y);
  }
  // Float Multiplication Operator Overload
  Vector2 operator*(const float& other) const
  {
    return Vector2(this->X *other, this->Y * other);
  }

  // Positional Variables
  float X;
  float Y;
};

// Structure: Vector3
//
// Description: A 3D Vector that Holds Positional Data
struct Vector3
{
  // Default Constructor
  Vector3()
  {
    X = 0.0f;
    Y = 0.0f;
    Z = 0.0f;
  }
  // Variable Set Constructor
  Vector3(float X_, float Y_, float Z_)
  {
    X = X_;
    Y = Y_;
    Z = Z_;
  }
  // Bool Equals Operator Overload
  bool operator==(const Vector3& other) const
  {
    return (this->X == other.X && this->Y == other.Y && this->Z == other.Z);
  }
  // Bool Not Equals Operator Overload
  bool operator!=(const Vector3& other) const
  {
    return !(this->X == other.X && this->Y == other.Y && this->Z == other.Z);
  }
  // Addition Operator Overload
  Vector3 operator+(const Vector3& right) const
  {
    return Vector3(this->X + right.X, this->Y + right.Y, this->Z + right.Z);
  }
  // Subtraction Operator Overload
  Vector3 operator-(const Vector3& right) const
  {
    return Vector3(this->X - right.X, this->Y - right.Y, this->Z - right.Z);
  }
  // Float Multiplication Operator Overload
  Vector3 operator*(const float& other) const
  {
    return Vector3(this->X * other, this->Y * other, this->Z * other);
  }
  // Float Division Operator Overload
  Vector3 operator/(const float& other) const
  {
    return Vector3(this->X / other, this->Y / other, this->Z / other);
  }

  // Positional Variables
  float X;
  float Y;
  float Z;
};

// Structure: Vertex
//
// Description: Model Vertex object that holds
//	a Position, Normal, and Texture Coordinate
struct Vertex
{
  // Position Vector
  Vector3 Position;

  // Normal Vector
  Vector3 Normal;

  // Texture Coordinate Vector
  Vector2 TextureCoordinate;
};

struct Material
{
  // Material Name
  std::string name;
  // Ambient Color
  Vector3 Ka;
  // Diffuse Color
  Vector3 Kd;
  // Specular Color
  Vector3 Ks;
  // Specular Exponent
  float Ns{30.0f};
  // Optical Density
  float Ni{0.0f};
  // Dissolve
  float d{1.0f};
  // Illumination
  int illum{0};
  // Ambient Texture Map
  std::string map_Ka;
  // Diffuse Texture Map
  std::string map_Kd;
  // Specular Texture Map
  std::string map_Ks;
  // Specular Hightlight Map
  std::string map_Ns;
  // Alpha Texture Map
  std::string map_d;
  // Bump Map
  std::string map_bump;
};

// Structure: Mesh
//
// Description: A Simple Mesh Object that holds
//	a name, a vertex list, and an index list
struct Mesh
{
  // Default Constructor
  Mesh()
  {

  }
  // Variable Set Constructor
  Mesh(std::vector<Vertex>& _Vertices, std::vector<unsigned int>& _Indices)
  {
    Vertices = _Vertices;
    Indices = _Indices;
  }
  // Mesh Name
  std::string MeshName;
  // Vertex List
  std::vector<Vertex> Vertices;
  // Index List
  std::vector<unsigned int> Indices;

  // Material
  Material MeshMaterial;
};

// Namespace: Math
//
// Description: The namespace that holds all of the math
//	functions need for OBJL
namespace math
{
// Vector3 Cross Product
Vector3 CrossV3(const Vector3 a, const Vector3 b)
{
  return Vector3(a.Y * b.Z - a.Z * b.Y,
                 a.Z * b.X - a.X * b.Z,
                 a.X * b.Y - a.Y * b.X);
}

// Vector3 Magnitude Calculation
float MagnitudeV3(const Vector3 in)
{
  return (sqrtf(powf(in.X, 2) + powf(in.Y, 2) + powf(in.Z, 2)));
}

// Vector3 DotProduct
float DotV3(const Vector3 a, const Vector3 b)
{
  return (a.X * b.X) + (a.Y * b.Y) + (a.Z * b.Z);
}

// Angle between 2 Vector3 Objects
float AngleBetweenV3(const Vector3 a, const Vector3 b)
{
  float angle = DotV3(a, b);
  angle /= (MagnitudeV3(a) * MagnitudeV3(b));
  return angle = acosf(angle);
}

// Projection Calculation of a onto b
Vector3 ProjV3(const Vector3 a, const Vector3 b)
{
  Vector3 bn = b / MagnitudeV3(b);
  return bn * DotV3(a, bn);
}
}

// Namespace: Algorithm
//
// Description: The namespace that holds all of the
// Algorithms needed for OBJL
namespace algorithm
{
// Vector3 Multiplication Opertor Overload
Vector3 operator*(const float& left, const Vector3& right)
{
  return Vector3(right.X * left, right.Y * left, right.Z * left);
}

// A test to see if P1 is on the same side as P2 of a line segment ab
bool SameSide(Vector3 p1, Vector3 p2, Vector3 a, Vector3 b)
{
  Vector3 cp1 = math::CrossV3(b - a, p1 - a);
  Vector3 cp2 = math::CrossV3(b - a, p2 - a);

  if (math::DotV3(cp1, cp2) >= 0)
    return true;
  else
    return false;
}

// Generate a cross produect normal for a triangle
Vector3 GenTriNormal(Vector3 t1, Vector3 t2, Vector3 t3)
{
  Vector3 u = t2 - t1;
  Vector3 v = t3 - t1;

  Vector3 normal = math::CrossV3(u,v);

  return normal;
}

// Check to see if a Vector3 Point is within a 3 Vector3 Triangle
bool inTriangle(Vector3 point, Vector3 tri1, Vector3 tri2, Vector3 tri3)
{
  // Test to see if it is within an infinite prism that the triangle outlines.
  bool within_tri_prisim = SameSide(point, tri1, tri2, tri3) && SameSide(point, tri2, tri1, tri3)
      && SameSide(point, tri3, tri1, tri2);

  // If it isn't it will never be on the triangle
  if (!within_tri_prisim)
    return false;

  // Calulate Triangle's Normal
  Vector3 n = GenTriNormal(tri1, tri2, tri3);

  // Project the point onto this normal
  Vector3 proj = math::ProjV3(point, n);

  // If the distance from the triangle to the point is 0
  //	it lies on the triangle
  if (math::MagnitudeV3(proj) == 0)
    return true;
  else
    return false;
}

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t\n\r")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

// Split a String into a string array at a given token
inline void split(const std::string &in,
                  std::vector<std::string> &out,
                  std::string token)
{
  out.clear();

  std::string temp;

  for (int i = 0; i < int(in.size()); i++)
  {
    std::string test = in.substr(i, token.size());

    if (test == token)
    {
      if (!temp.empty())
      {
        out.push_back(temp);
        temp.clear();
        i += int(token.size()) - 1;
      }
      else
      {
        out.push_back("");
      }
    }
    else if (i + token.size() >= in.size())
    {
      temp += in.substr(i, token.size());
      out.push_back(temp);
      break;
    }
    else
    {
      temp += in[i];
    }
  }
}

// Get tail of string after first token and possibly following spaces
inline std::string tail(const std::string &in)
{
  size_t token_start = in.find_first_not_of(" \t");
  size_t space_start = in.find_first_of(" \t", token_start);
  size_t tail_start = in.find_first_not_of(" \t", space_start);
  size_t tail_end = in.find_last_not_of(" \t");
  if (tail_start != std::string::npos && tail_end != std::string::npos)
  {
    return in.substr(tail_start, tail_end - tail_start + 1);
  }
  else if (tail_start != std::string::npos)
  {
    return in.substr(tail_start);
  }
  return "";
}

// Get first token of string
inline std::string firstToken(const std::string &in)
{
  if (!in.empty())
  {
    size_t token_start = in.find_first_not_of(" \t");
    size_t token_end = in.find_first_of(" \t", token_start);
    if (token_start != std::string::npos && token_end != std::string::npos)
    {
      return in.substr(token_start, token_end - token_start);
    }
    else if (token_start != std::string::npos)
    {
      return in.substr(token_start);
    }
  }
  return "";
}

// Get element at given index position
template <class T>
inline const T & getElement(const std::vector<T> &elements, std::string &index)
{
  int idx = std::stoi(index);
  if (idx < 0)
    idx = int(elements.size()) + idx;
  else
    idx--;
  return elements[idx];
}
}

// Class: Loader
//
// Description: The OBJ Model Loader
class Loader
{
public:

  //################################################################################################
  // Load a file into the loader
  //
  // If file is loaded return true
  //
  // If the file is unable to be found
  // or unable to be loaded return false
  bool LoadFile(const std::string& path)
  {
    // If the file is not an .obj file return false
    if (path.substr(path.size() - 4, 4) != ".obj")
      return false;

    std::ifstream objStream(path);

    if (!objStream.is_open())
      return false;

    return LoadStreams(objStream, [&](const std::string& fileName, const std::function<void(std::istream&)>& closure)
    {
      // Generate a path to the material file
      std::vector<std::string> temp;
      algorithm::split(path, temp, "/");

      std::string pathtomat = "";

      if (temp.size() != 1)
        for (size_t i = 0; i < temp.size() - 1; i++)
          pathtomat += temp[i] + "/";

      pathtomat += fileName;

#ifdef OBJL_CONSOLE_OUTPUT
      std::cout << "Load material: '" << pathtomat << "'" << std::endl;
#endif

      std::cout << "pathtomat: " << pathtomat << std::endl;
      std::ifstream mtlStream(pathtomat);

      if(!mtlStream.is_open())
        return;

      closure(mtlStream);
    });
  }

  //################################################################################################
  bool LoadStreams(std::istream& objStream, const std::function<void(const std::string&, const std::function<void(std::istream&)>&)>& loadMTL)
  {

    LoadedMeshes.clear();
    LoadedVertices.clear();
    LoadedIndices.clear();

    std::vector<Vector3> Positions;
    std::vector<Vector2> TCoords;
    std::vector<Vector3> Normals;

    std::vector<Vertex> Vertices;
    std::vector<unsigned int> Indices;
    std::string matName;

    //std::vector<std::string> MeshMatNames;

    bool listening = false;
    std::string meshname;

    Mesh tempMesh;

#ifdef OBJL_CONSOLE_OUTPUT
    const unsigned int outputEveryNth = 1000;
    unsigned int outputIndicator = outputEveryNth;
#endif

    std::string curline;
    while (std::getline(objStream, curline))
    {
      curline = algorithm::trim(curline);
      std::string firstToken = algorithm::firstToken(curline);


#ifdef OBJL_CONSOLE_OUTPUT
      if ((outputIndicator = ((outputIndicator + 1) % outputEveryNth)) == 1)
      {
        if (!meshname.empty())
        {
          std::cout
              << "\r- " << meshname
              << "\t| vertices > " << Positions.size()
              << "\t| texcoords > " << TCoords.size()
              << "\t| normals > " << Normals.size()
              << "\t| triangles > " << (Vertices.size() / 3);
        }
      }
#endif

      // Generate a Mesh Object or Prepare for an object to be created
      if(firstToken == "o" || firstToken == "g") // || curline[0] == 'g'
      {
        if (!listening)
        {
          listening = true;

          //if (firstToken == "o" || firstToken == "g")
          //{
            meshname = algorithm::tail(curline);
          //}
          //else
          //{
          //  meshname = "unnamed";
          //}
        }
        else
        {
          // Generate the mesh to put into the array

          if (!Indices.empty() && !Vertices.empty())
          {
            // Create Mesh
            tempMesh = Mesh(Vertices, Indices);
            tempMesh.MeshName = meshname;
            tempMesh.MeshMaterial.name = matName;

            // Insert Mesh
            LoadedMeshes.push_back(tempMesh);

            // Cleanup
            Vertices.clear();
            Indices.clear();
            meshname.clear();
            matName.clear();

            meshname = algorithm::tail(curline);
          }
          else
          {
            //if (firstToken == "o" || firstToken == "g")
            //{
              meshname = algorithm::tail(curline);
            //}
            //else
            //{
            //  meshname = "unnamed";
            //}
          }
        }
#ifdef OBJL_CONSOLE_OUTPUT
        std::cout << std::endl;
        outputIndicator = 0;
#endif
      }

      // Generate a Vertex Position
      if (firstToken == "v")
      {
        std::vector<std::string> spos;
        Vector3 vpos;
        algorithm::split(algorithm::tail(curline), spos, " ");

        vpos.X = readFloat(spos[0]);
        vpos.Y = readFloat(spos[1]);
        vpos.Z = readFloat(spos[2]);

        Positions.push_back(vpos);
      }
      // Generate a Vertex Texture Coordinate
      if (firstToken == "vt")
      {
        std::vector<std::string> stex;
        Vector2 vtex;
        algorithm::split(algorithm::tail(curline), stex, " ");

        vtex.X = readFloat(stex[0]);
        vtex.Y = readFloat(stex[1]);

        TCoords.push_back(vtex);
      }
      // Generate a Vertex Normal;
      if (firstToken == "vn")
      {
        std::vector<std::string> snor;
        Vector3 vnor;
        algorithm::split(algorithm::tail(curline), snor, " ");

        vnor.X = readFloat(snor[0]);
        vnor.Y = readFloat(snor[1]);
        vnor.Z = readFloat(snor[2]);

        Normals.push_back(vnor);
      }
      // Generate a Face (vertices & indices)
      if (firstToken == "f")
      {
        // Generate the vertices
        std::vector<Vertex> vVerts;
        GenVerticesFromRawOBJ(vVerts, Positions, TCoords, Normals, curline);

        // Add Vertices
        for (int i = 0; i < int(vVerts.size()); i++)
        {
          Vertices.push_back(vVerts[i]);

          LoadedVertices.push_back(vVerts[i]);
        }

        std::vector<unsigned int> iIndices;

        VertexTriangluation(iIndices, vVerts);

        // Add Indices
        for (int i = 0; i < int(iIndices.size()); i++)
        {
          unsigned int indnum = unsigned((Vertices.size()) - vVerts.size()) + iIndices[i];
          Indices.push_back(indnum);

          indnum = unsigned((LoadedVertices.size()) - vVerts.size()) + iIndices[i];
          LoadedIndices.push_back(indnum);

        }
      }

      // Get Mesh Material Name
      if(firstToken == "usemtl")
      {
        matName = algorithm::tail(curline);
        std::cout << "usemtl: '" << matName << "' size: " << matName.size() << std::endl;
        //MeshMatNames.push_back(algorithm::tail(curline));

        // Create new Mesh, if Material changes within a group
        if (!Indices.empty() && !Vertices.empty())
        {
          // Create Mesh
          tempMesh = Mesh(Vertices, Indices);
          tempMesh.MeshName = meshname;
          tempMesh.MeshMaterial.name = matName;

          int i = 2;
          while(1) {
            tempMesh.MeshName = meshname + "_" + std::to_string(i);

            for (auto &m : LoadedMeshes)
              if (m.MeshName == tempMesh.MeshName)
                continue;
            break;
          }

          // Insert Mesh
          LoadedMeshes.push_back(tempMesh);

          // Cleanup
          Vertices.clear();
          Indices.clear();
          matName.clear();
        }

#ifdef OBJL_CONSOLE_OUTPUT
        outputIndicator = 0;
#endif
      }

      // Load Materials
      if (firstToken == "mtllib")
      {        
        std::cout << "mtllib" << std::endl;
        loadMTL(algorithm::tail(curline), [&](std::istream& mtlStream)
        {
          std::cout << "loadMTL " << std::endl;
          LoadMaterials(mtlStream);
        });
      }
    }

#ifdef OBJL_CONSOLE_OUTPUT
    std::cout << std::endl;
#endif

    // Deal with last mesh

    if (!Indices.empty() && !Vertices.empty())
    {
      // Create Mesh
      tempMesh = Mesh(Vertices, Indices);
      tempMesh.MeshName = meshname;
      tempMesh.MeshMaterial.name = matName;

      // Insert Mesh
      LoadedMeshes.push_back(tempMesh);
    }

    // Set Materials for each Mesh
    for (size_t i = 0; i < LoadedMeshes.size(); i++)
    {
      //std::string matname = MeshMatNames[i];
      const std::string& matname = LoadedMeshes[i].MeshMaterial.name;
      std::cout << "Find material: " << matname << std::endl;

      // Find corresponding material name in loaded materials
      // when found copy material variables into mesh material
      for (size_t j = 0; j < LoadedMaterials.size(); j++)
      {
        std::cout << "    - " << LoadedMaterials[j].name << std::endl;
        if (LoadedMaterials[j].name == matname)
        {
          LoadedMeshes[i].MeshMaterial = LoadedMaterials[j];
          break;
        }
      }
    }

    if (LoadedMeshes.empty() && LoadedVertices.empty() && LoadedIndices.empty())
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  // Loaded Mesh Objects
  std::vector<Mesh> LoadedMeshes;
  // Loaded Vertex Objects
  std::vector<Vertex> LoadedVertices;
  // Loaded Index Positions
  std::vector<unsigned int> LoadedIndices;
  // Loaded Material Objects
  std::vector<Material> LoadedMaterials;

private:
  // Generate vertices from a list of positions,
  //	tcoords, normals and a face line
  void GenVerticesFromRawOBJ(std::vector<Vertex>& oVerts,
                             const std::vector<Vector3>& iPositions,
                             const std::vector<Vector2>& iTCoords,
                             const std::vector<Vector3>& iNormals,
                             std::string icurline)
  {
    std::vector<std::string> sface, svert;
    Vertex vVert;
    algorithm::split(algorithm::tail(icurline), sface, " ");

    bool noNormal = false;

    // For every given vertex do this
    for (int i = 0; i < int(sface.size()); i++)
    {
      // To handle this situation where we have a trailing space "f 1/1/1 2/2/1 3/3/2 4/4/2 \r"
      if(sface[i] == "\r")
        break;

      // See What type the vertex is.
      int vtype{0};

      algorithm::split(sface[i], svert, "/");

      // Check for just position - v1
      if (svert.size() == 1)
      {
        // Only position
        vtype = 1;
      }

      // Check for position & texture - v1/vt1
      if (svert.size() == 2)
      {
        // Position & Texture
        vtype = 2;
      }

      // Check for Position, Texture and Normal - v1/vt1/vn1
      // or if Position and Normal - v1//vn1
      if (svert.size() == 3)
      {
        if (svert[1] != "")
        {
          // Position, Texture, and Normal
          vtype = 4;
        }
        else
        {
          // Position & Normal
          vtype = 3;
        }
      }

      // Calculate and store the vertex
      try
      {
        switch (vtype)
        {
        case 1: // P
        {
          vVert.Position = algorithm::getElement(iPositions, svert[0]);
          vVert.TextureCoordinate = Vector2(0, 0);
          noNormal = true;
          oVerts.push_back(vVert);
          break;
        }
        case 2: // P/T
        {
          vVert.Position = algorithm::getElement(iPositions, svert[0]);
          vVert.TextureCoordinate = algorithm::getElement(iTCoords, svert[1]);
          noNormal = true;
          oVerts.push_back(vVert);
          break;
        }
        case 3: // P//N
        {
          vVert.Position = algorithm::getElement(iPositions, svert[0]);
          vVert.TextureCoordinate = Vector2(0, 0);
          vVert.Normal = algorithm::getElement(iNormals, svert[2]);
          oVerts.push_back(vVert);
          break;
        }
        case 4: // P/T/N
        {
          vVert.Position = algorithm::getElement(iPositions, svert[0]);
          vVert.TextureCoordinate = algorithm::getElement(iTCoords, svert[1]);
          vVert.Normal = algorithm::getElement(iNormals, svert[2]);
          oVerts.push_back(vVert);
          break;
        }
        default:
        {
          break;
        }
        }
      } catch (...)
      {
        std::cerr << "Error loading OBJ file at " << __FILE__ << ":" << __LINE__ << std::endl;
      }
    }

    // take care of missing normals
    // these may not be truly acurate but it is the
    // best they get for not compiling a mesh with normals
    if (noNormal)
    {
      Vector3 A = oVerts[0].Position - oVerts[1].Position;
      Vector3 B = oVerts[2].Position - oVerts[1].Position;

      Vector3 normal = math::CrossV3(A, B);

      for (int i = 0; i < int(oVerts.size()); i++)
      {
        oVerts[i].Normal = normal;
      }
    }
  }

  // Triangulate a list of vertices into a face by printing
  //	inducies corresponding with triangles within it
  void VertexTriangluation(std::vector<unsigned int>& oIndices,
                           const std::vector<Vertex>& iVerts)
  {
    // If there are 2 or less verts,
    // no triangle can be created,
    // so exit
    if (iVerts.size() < 3)
    {
      return;
    }
    // If it is a triangle no need to calculate it
    if (iVerts.size() == 3)
    {
      oIndices.push_back(0);
      oIndices.push_back(1);
      oIndices.push_back(2);
      return;
    }

    // If it is a triangle no need to calculate it
    if (iVerts.size() == 4)
    {
      oIndices.push_back(0);
      oIndices.push_back(1);
      oIndices.push_back(2);
      oIndices.push_back(0);
      oIndices.push_back(2);
      oIndices.push_back(3);
      return;
    }

    // Create a list of vertices
    std::vector<Vertex> tVerts = iVerts;

    while (true)
    {
      // For every vertex
      for (size_t i = 0; i < tVerts.size(); i++)
      {
        // pPrev = the previous vertex in the list
        Vertex pPrev;
        if (i == 0)
        {
          pPrev = tVerts[tVerts.size() - 1];
        }
        else
        {
          pPrev = tVerts[i - 1];
        }

        // pCur = the current vertex;
        Vertex pCur = tVerts[i];

        // pNext = the next vertex in the list
        Vertex pNext;
        if (i == tVerts.size() - 1)
        {
          pNext = tVerts[0];
        }
        else
        {
          pNext = tVerts[i + 1];
        }

        // Check to see if there are only 3 verts left
        // if so this is the last triangle
        if (tVerts.size() == 3)
        {
          // Create a triangle from pCur, pPrev, pNext
          for (int j = 0; j < int(tVerts.size()); j++)
          {
            if (iVerts[j].Position == pCur.Position)
              oIndices.push_back(j);
            if (iVerts[j].Position == pPrev.Position)
              oIndices.push_back(j);
            if (iVerts[j].Position == pNext.Position)
              oIndices.push_back(j);
          }

          tVerts.clear();
          break;
        }
        if (tVerts.size() == 4)
        {
          // Create a triangle from pCur, pPrev, pNext
          for (int j = 0; j < int(iVerts.size()); j++)
          {
            if (iVerts[j].Position == pCur.Position)
              oIndices.push_back(j);
            if (iVerts[j].Position == pPrev.Position)
              oIndices.push_back(j);
            if (iVerts[j].Position == pNext.Position)
              oIndices.push_back(j);
          }

          Vector3 tempVec;
          for (int j = 0; j < int(tVerts.size()); j++)
          {
            if (tVerts[j].Position != pCur.Position
                && tVerts[j].Position != pPrev.Position
                && tVerts[j].Position != pNext.Position)
            {
              tempVec = tVerts[j].Position;
              break;
            }
          }

          // Create a triangle from pCur, pPrev, pNext
          for (int j = 0; j < int(iVerts.size()); j++)
          {
            if (iVerts[j].Position == pPrev.Position)
              oIndices.push_back(j);
            if (iVerts[j].Position == pNext.Position)
              oIndices.push_back(j);
            if (iVerts[j].Position == tempVec)
              oIndices.push_back(j);
          }

          tVerts.clear();
          break;
        }

        continue;

        // If Vertex is not an interior vertex
        float angle = math::AngleBetweenV3(pPrev.Position - pCur.Position, pNext.Position - pCur.Position) * (180.0f / 3.14159265359f);
        if (angle <= 0 && angle >= 180)
          continue;

        // If any vertices are within this triangle
        bool inTri = false;
        for (int j = 0; j < int(iVerts.size()); j++)
        {
          if (algorithm::inTriangle(iVerts[j].Position, pPrev.Position, pCur.Position, pNext.Position)
              && iVerts[j].Position != pPrev.Position
              && iVerts[j].Position != pCur.Position
              && iVerts[j].Position != pNext.Position)
          {
            inTri = true;
            break;
          }
        }
        if (inTri)
          continue;

        // Create a triangle from pCur, pPrev, pNext
        for (int j = 0; j < int(iVerts.size()); j++)
        {
          if (iVerts[j].Position == pCur.Position)
            oIndices.push_back(j);
          if (iVerts[j].Position == pPrev.Position)
            oIndices.push_back(j);
          if (iVerts[j].Position == pNext.Position)
            oIndices.push_back(j);
        }

        // Delete pCur from the list
        for (int j = 0; j < int(tVerts.size()); j++)
        {
          if (tVerts[j].Position == pCur.Position)
          {
            tVerts.erase(tVerts.begin() + j);
            break;
          }
        }

        // reset i to the start
        // -1 since loop will add 1 to it
        i = -1;
      }

      // if no triangles were created
      if (oIndices.size() == 0)
        break;

      // if no more vertices
      if (tVerts.size() == 0)
        break;
    }
  }

  // Load Materials from .mtl file
  bool LoadMaterials(std::istream& inputStream)
  {
    Material tempMaterial;

    bool listening = false;

    // Go through each line looking for material variables
    std::string curline;
    while (std::getline(inputStream, curline))
    {
      auto a = curline;
      curline = algorithm::trim(curline);

      std::cout << "A: \"" << a << "\" B: \"" << curline << "\"" << std::endl;

      std::string firstToken = algorithm::firstToken(curline);

      // new material and material name
      if (firstToken == "newmtl")
      {
        if (!listening)
        {
          listening = true;

          if (curline.size() > 7)
          {            
            tempMaterial.name = algorithm::tail(curline);            

          }
          else
          {
            tempMaterial.name = "none";
          }
        }
        else
        {
          // Generate the material

          // Push Back loaded Material
          LoadedMaterials.push_back(tempMaterial);

          // Clear Loaded Material
          tempMaterial = Material();

          if (curline.size() > 7)
          {
            tempMaterial.name = algorithm::tail(curline);
          }
          else
          {
            tempMaterial.name = "none";
          }          
        }

        std::cout << "newmtl: '" << tempMaterial.name << "' size: " << tempMaterial.name.size() << std::endl;
      }
      // Ambient Color
      if (firstToken == "Ka")
      {
        std::vector<std::string> temp;
        algorithm::split(algorithm::tail(curline), temp, " ");

        if (temp.size() != 3)
          continue;

        tempMaterial.Ka.X = readFloat(temp[0]);
        tempMaterial.Ka.Y = readFloat(temp[1]);
        tempMaterial.Ka.Z = readFloat(temp[2]);
      }
      // Diffuse Color
      if (firstToken == "Kd")
      {
        std::vector<std::string> temp;
        algorithm::split(algorithm::tail(curline), temp, " ");

        if (temp.size() != 3)
          continue;

        tempMaterial.Kd.X = readFloat(temp[0]);
        tempMaterial.Kd.Y = readFloat(temp[1]);
        tempMaterial.Kd.Z = readFloat(temp[2]);
      }
      // Specular Color
      if (firstToken == "Ks")
      {
        std::vector<std::string> temp;
        algorithm::split(algorithm::tail(curline), temp, " ");

        if (temp.size() != 3)
          continue;

        tempMaterial.Ks.X = readFloat(temp[0]);
        tempMaterial.Ks.Y = readFloat(temp[1]);
        tempMaterial.Ks.Z = readFloat(temp[2]);
      }
      // Specular Exponent
      if (firstToken == "Ns")
      {
        tempMaterial.Ns = readFloat(algorithm::tail(curline));
      }
      // Optical Density
      if (firstToken == "Ni")
      {
        tempMaterial.Ni = readFloat(algorithm::tail(curline));
      }
      // Dissolve
      if (firstToken == "d")
      {
        tempMaterial.d = readFloat(algorithm::tail(curline));
      }
      // Illumination
      if (firstToken == "illum")
      {
        tempMaterial.illum = std::stoi(algorithm::tail(curline));
      }
      // Ambient Texture Map
      if (firstToken == "map_Ka")
      {
        tempMaterial.map_Ka = algorithm::tail(curline);
      }
      // Diffuse Texture Map
      if (firstToken == "map_Kd")
      {
        tempMaterial.map_Kd = algorithm::tail(curline);
      }
      // Specular Texture Map
      if (firstToken == "map_Ks")
      {
        tempMaterial.map_Ks = algorithm::tail(curline);
      }
      // Specular Hightlight Map
      if (firstToken == "map_Ns")
      {
        tempMaterial.map_Ns = algorithm::tail(curline);
      }
      // Alpha Texture Map
      if (firstToken == "map_d")
      {
        tempMaterial.map_d = algorithm::tail(curline);
      }
      // Bump Map
      if (firstToken == "map_Bump" || firstToken == "map_bump" || firstToken == "bump")
      {
        tempMaterial.map_bump = algorithm::tail(curline);
      }
    }

    // Deal with last material

    // Push Back loaded Material
    LoadedMaterials.push_back(tempMaterial);

    // Test to see if anything was loaded
    // If not return false
    if (LoadedMaterials.empty())
      return false;
    // If so return true
    else
      return true;
  }
};
}
