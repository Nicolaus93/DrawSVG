#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = svg_2_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate; 
  size_t supersample_size = 4 * (this->target_w * this->sample_rate) * (this->target_h * this->sample_rate);
  this->sample_buffer.resize(supersample_size);
  fill(this->sample_buffer.begin(), this->sample_buffer.end(), 255);
}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;
  
  // Task 4: 
  // You may want to modify this for supersampling support
  size_t supersample_size = 4 * this->target_w * this->target_h;
  this->sample_buffer = std::vector<unsigned char>(supersample_size, 255);
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= target_w ) return;
  if ( sy < 0 || sy >= target_h ) return;

  // fill sample - NOT doing alpha blending!
  render_target[4 * (sx + sy * target_w)    ] = (uint8_t) (color.r * 255);
  render_target[4 * (sx + sy * target_w) + 1] = (uint8_t) (color.g * 255);
  render_target[4 * (sx + sy * target_w) + 2] = (uint8_t) (color.b * 255);
  render_target[4 * (sx + sy * target_w) + 3] = (uint8_t) (color.a * 255);

}


void SoftwareRendererImp::update_sample_buffer(float x, float y, Color color) {
  int sx = (int) floor(x * this->sample_rate);
  int sy = (int) floor(y * this->sample_rate);

  if ( sx < 0 || sx >= target_w * this->sample_rate) return;
  if ( sy < 0 || sy >= target_h * this->sample_rate) return;

  sample_buffer[4 * (sx + sy * target_w * this->sample_rate)    ] = (uint8_t) (color.r * 255);
  sample_buffer[4 * (sx + sy * target_w * this->sample_rate) + 1] = (uint8_t) (color.g * 255);
  sample_buffer[4 * (sx + sy * target_w * this->sample_rate) + 2] = (uint8_t) (color.b * 255);
  sample_buffer[4 * (sx + sy * target_w * this->sample_rate) + 3] = (uint8_t) (color.a * 255);
}


void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 2: 
  // Implement line rasterization
  float dx = (x1 - x0);
  float dy = (y1 - y0);
  float m = dx != 0? dy / dx : 1;
  int eps = 0;

  if (dy == 0) {
    for (int x = (int) std::min(x0, x1); x < std::max(x0, x1); ++x ) {
      this->rasterize_point(x, y0, color);
    }
    return;
  }

  if ( m > 0 && m < 1 ) {
    if (dx < 0) {
      this->rasterize_line(x1, y1, x0, y0, color);
      return;
    }
    float y{y0};
    for (int x = (int) x0; x <= (int) x1; ++x) {
      this->rasterize_point(x, y, color);
      eps += dy;
      if ( (eps << 1) >= dx ) {
        ++y;
        eps -= dx;
      }
    }
  }

  else if (m >= 1) {
    if (dy < 0) {
      this->rasterize_line(x1, y1, x0, y0, color);
      return;
    }
    int x{x0};
    for (int y = (int) y0; y <= (int) y1; ++y) {
      this->rasterize_point(x, y, color);
      eps += dx;
      if ( (eps << 1) >= dy ) {
        ++x;
        eps -= dy;
      }
    }
  }

  else if (m < 0 && m > -1) {
    if (dx < 0) {
      this->rasterize_line(x1, y1, x0, y0, color);
      return;
    }
    int y{y0};
    for (int x = (int) x0; x <= x1; ++x) {
      this->rasterize_point(x, y, color);
      eps += dy;
      if ( (eps << 1) < -dx ) {
        --y;
        eps += dx;
      }
    }
  }

  else {
    if ( dy < 0 ) {
      this->rasterize_line(x1, y1, x0, y0, color);
      return;
    }
    int x{x0};
    for (int y = (int) y0; y < y1; ++y) {
      this->rasterize_point(x, y, color);
      eps += dx;
      if ( (eps << 1 ) < -dy ) {
        --x;
        eps += dy;
      }
    }
  }
}

/*
checks if a point is in a triangle
https://stackoverflow.com/questions/1585459/whats-the-most-efficient-way-to-detect-triangle-triangle-intersections
*/ 
bool is_pnt_in_tri(const std::vector<Vector2D>& triangle, const Vector2D& p) {
  Vector2D ab = triangle.at(1) - triangle.at(0);
  Vector2D bc = triangle.at(2) - triangle.at(1);
  Vector2D ca = triangle.at(0) - triangle.at(2);
  Vector2D ap = p - triangle.at(0);
  Vector2D bp = p - triangle.at(1);
  Vector2D cp = p - triangle.at(2);
  if (cross(ab, ap) > 0 && cross(bc, bp) > 0 && cross(ca, cp) > 0 ) {
    return true;
  }
  return false;
}


// check http://thirdpartyninjas.com/blog/2008/10/07/line-segment-intersection/
bool line_intersection(Vector2D p1, Vector2D p2, Vector2D p3, Vector2D p4) {
  float den = (p4.y - p3.y) * (p2.x - p1.x) - (p4.x - p3.x) * (p2.y - p1.y);
  if (den == 0) {
    return false;
  }
  float t = ((p4.x - p3.x) * (p1.y - p3.y) - (p4.y - p3.y) * (p1.x - p3.x)) / den;
  if (t > 1 || t < 0) {
    return false;
  }
  float s = ((p2.x - p1.x) * (p1.y - p3.y) - (p2.y - p1.y) * (p1.x - p3.x)) / den;
  if (s > 1 || s < 0) {
    return false;
  }
  return true;
}


bool triangle_int(const std::vector<Vector2D>& tri1, 
                  const std::vector<Vector2D>& tri2) {
  // check for vertex inside other tri
  for (auto& p : tri1) {
    if (is_pnt_in_tri(tri2, p)) {
      return true;
    }
  }
  for (auto& p : tri2) {
    if (is_pnt_in_tri(tri1, p)) {
      return true;
    }
  }

  // check for sides intersections
  std::vector<std::pair<Vector2D, Vector2D>> sides_1{
    std::make_pair(tri1.at(0), tri1.at(1)), 
    std::make_pair(tri1.at(1), tri1.at(2)), 
    std::make_pair(tri1.at(2), tri1.at(0))
  };
  std::vector<std::pair<Vector2D, Vector2D>> sides_2{
    std::make_pair(tri2.at(0), tri2.at(1)), 
    std::make_pair(tri2.at(1), tri2.at(2)), 
    std::make_pair(tri2.at(2), tri2.at(0))
  };
  for (auto& side1 : sides_1) {
    for (auto& side2 : sides_2) {
      if (line_intersection(side1.first, side1.second, side2.first, side2.second)) {
        return true;
      }
    }
  }

  return false;
}


bool box_triangle_int(const Rect& bbox, 
                      const std::vector<Vector2D>& triangle) {
  float min_x = bbox.position.x;
  float max_x = min_x + bbox.dimension.x;
  float min_y = bbox.position.y;
  float max_y = min_y + bbox.dimension.y;

  std::vector<Vector2D> upper_tri{Vector2D(min_x, min_y), Vector2D(max_x, min_y), Vector2D(min_x, max_y)};
  std::vector<Vector2D> lower_tri{Vector2D(max_x, max_y), Vector2D(max_x, min_y), Vector2D(min_x, max_y)};
  if (!triangle_int(upper_tri, triangle) && !triangle_int(lower_tri, triangle)) {
    return false;
  }
  return true;
}


std::vector<Vector2D> get_neighbours(const Vector2D& pnt, int sample_rate) {
  std::vector<Vector2D> neighbours;
  for (int i = 0; i < sample_rate; ++i) {
    for (int j = 0; j < sample_rate; ++j) {
      neighbours.emplace_back(Vector2D(pnt.x + (float) i / sample_rate, pnt.y + (float) j / sample_rate));
    }
  }
  return neighbours;
}


/*recursive function -> use level for debugging*/
void find_raster_points( std::vector<Vector2D>& triangle,
                         Rect& bbox,
                         std::vector<Vector2D>& raster_pts,
                         int sample_rate,
                         int level = 0) {
  
  float w = bbox.dimension.x;
  float h = bbox.dimension.y;
  if (w > 10 || h > 10) {
    if (!box_triangle_int(bbox, triangle)) {
      return;
    }
    // if still big, split box
    std::vector<Vector2D> positions{
      Vector2D(), 
      Vector2D{bbox.dimension.x / 2, 0}, 
      bbox.dimension / 2,
      Vector2D{0, bbox.dimension.y / 2},
    };
    Rect new_box = Rect();
    new_box.dimension = bbox.dimension / 2;
    for (auto& offset : positions) {
      new_box.position = bbox.position + offset;
      find_raster_points(triangle, new_box, raster_pts, sample_rate, level + 1);  
    }

  } else {

    if (!box_triangle_int(bbox, triangle)) {
      return;
    }

    float min_x = bbox.position.x;
    float max_x = min_x + bbox.dimension.x;
    float min_y = bbox.position.y;
    float max_y = min_y + bbox.dimension.y;
    for (int x{min_x}; x < max_x; ++x) {
      for (int y{min_y}; y < max_y; ++y) {
        std::vector<Vector2D> pts_to_check = get_neighbours(Vector2D(x, y), sample_rate);
        for (auto& pnt : pts_to_check) {
          if (is_pnt_in_tri(triangle, pnt)) {
            raster_pts.emplace_back(pnt);
          }
        }
      }
    }
  }
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {

  std::vector<Vector2D> triangle{Vector2D(x0, y0), Vector2D(x1, y1), Vector2D(x2, y2)};
  float min_x = std::min(std::min(x0, x1), x2);
  float max_x = std::max(std::max(x0, x1), x2);
  float min_y = std::min(std::min(y0, y1), y2);
  float max_y = std::max(std::max(y0, y1), y2);
  Rect bbox = Rect();
  bbox.position = Vector2D(min_x, min_y);
  bbox.dimension = Vector2D(max_x - min_x, max_y - min_y);

  std::vector<Vector2D> sample_pts;
  find_raster_points(triangle, bbox, sample_pts, this->sample_rate);
  for (auto& p : sample_pts) {
    this->update_sample_buffer(p.x, p.y, color);
  }
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization

}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  for (int x = 0; x < this->target_w; ++x) {
    for (int y = 0; y < this->target_h; ++y) {
      std::vector<Vector2D> pts_to_check = get_neighbours(Vector2D(x, y), this->sample_rate);
      // iterate over rgba (0-4)
      for ( int i = 0; i < 4; ++i ) {
        int sum = 0;
        int count = 0;
        for (auto& p : pts_to_check) {
          int sx = (int) floor(p.x * this->sample_rate);
          int sy = (int) floor(p.y * this->sample_rate);

          if ( sx < 0 || sx >= target_w * this->sample_rate) continue;
          if ( sy < 0 || sy >= target_h * this->sample_rate) continue;

          count += 1;
          sum += sample_buffer.at(4 * (sx + sy * target_w * this->sample_rate) + i);
        }

        render_target[4 * (x + y * target_w) + i] = sum / count;
      }
    }
  }
  
}


} // namespace CMU462
