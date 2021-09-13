#include <nori/object.h>
#include <nori/texture.h>
#include <stb_image.h>

NORI_NAMESPACE_BEGIN
template <typename T>
class ImageTexture : public Texture<T> {
public:
    ImageTexture(const PropertyList& props);

    virtual std::string toString() const override;

    virtual T eval(const Point2f& uv) override {
        if (!m_load_texture_successed) {
            //std::cout << "The image texture did not successfully loaded!";
            return Color3f(1.0f, 0.0f, 0.0f);
        }
        float extend_u = (uv.x() - (int)uv.x()) * (m_width - 1), extend_v = (uv.y()-(int)uv.y()) * (m_height - 1);
        int u_left = (int)extend_u, u_right = (int)(extend_u + 1);
        int v_down = (int)extend_v, v_up = (int)(extend_v + 1);
        float weight_left = u_right - extend_u, weight_right = extend_u - u_left;
        float weight_down = v_up - extend_v, weight_up = extend_v - v_down;

        //Seek pixel index
        int idx_l_d = m_channels * (m_width * (v_down) + u_left);
        int idx_l_u = m_channels * (m_width * (v_up % m_height) + u_left);
        int idx_r_d = m_channels * (m_width * (v_down) + u_right % m_width);
        int idx_r_u = m_channels * (m_width * (v_up % m_height) + u_right % m_width);
        if (idx_r_u > m_channels * m_width * m_height) {
            std::cout << "index out of the range! " << std::endl;
            std::cout << "UV: " << uv.x() << ", " << uv.y() << std::endl;
            return Color3f(1.0f);
        }
        unsigned char ref_0 = 0;
        return weight_left * weight_down * Color3f(m_pixels[idx_l_d] / 255.0f, m_pixels[idx_l_d + 1] / 255.0f, m_pixels[idx_l_d + 2] / 255.0f) + \
            weight_left * weight_up * Color3f(m_pixels[idx_l_u] / 255.0f, m_pixels[idx_l_u + 1] /255.0f, m_pixels[idx_l_u + 2] / 255.0f) + \
            weight_right * weight_down * Color3f(m_pixels[idx_r_d] / 255.0f, m_pixels[idx_r_d + 1] / 255.0f, m_pixels[idx_r_d + 2] / 255.0f) + \
            weight_right * weight_up * Color3f(m_pixels[idx_r_u] / 255.0f, m_pixels[idx_r_u + 1] / 255.0f, m_pixels[idx_r_u + 2] / 255.0f);
    }

    /*
    bool Load_BMP_Texture(std::string filename, int size, int& width, int& height) // Creates Texture From A Bitmap File
    {
        if (size > 4 || size < 3) {
            std::cout << "The size for bmp file only support 3 or 4!";
            return false;
        }
        FILE* fp;
        int total_bytes;

        fp = fopen(filename.c_str(), "rb");
        if (fp == 0) {
            std::cout << "Can not open the image texture file! " << std::endl;
            return false;
        }
        fseek(fp, 0x0012, SEEK_SET); //定位指针位置到：SEEK_SET(0) + 0x0012(18)处
        fread(&width, 4, 1, fp); width = abs(width);//void *width指向一个带有最小尺寸为4(读取的每个元素为4字节)*1(一次读取1个元素)内存块的指针
        fread(&height, 4, 1, fp); height = abs(height);//这里的height读取出来是-512
        const int BMP_Header_Length = 54;
        fseek(fp, BMP_Header_Length, SEEK_SET);//定位指针位置到：SEEK_SET(0) + BMP_Header_Length(54)处

        int line_bytes = width * size; //这是一个RGBA类型的图像，因此需要进行位深度的绑定，也就是4个Byte，32个位。

        total_bytes = line_bytes * abs(height);//总的Byte数：width * height * 位深度

        m_pixels = (unsigned char*)malloc(total_bytes);//指针的存储单位是Byte
        if (!m_pixels) {
            cout << "[Image Texture] allocate space failed! " << std::endl;
            fclose(fp);
            return false;
        }
        int success_read_num = fread(m_pixels, 1, total_bytes, fp);//fread会返回成功读取的个数
        if (success_read_num <= 0) {
            cout << "[Image Texture] Can not read pixels! " << endl;
            free(m_pixels);
            fclose(fp);
            return false;
        }
        fclose(fp);
        std::cout << "Load Texture: " << filename << " successfully! ";
        std::cout << "With pixels: " << success_read_num;
        std::cout << ", With width: " << width;
        std::cout << ", With height: " << height << std::endl;
        return true;
    }
    */

    virtual ~ImageTexture() {
        free(m_pixels);
    }

protected:
    std::string m_filepath;
    unsigned char* m_pixels;
    int m_width, m_height, m_channels;
    bool m_load_texture_successed;
};

template <>
ImageTexture<Color3f>::ImageTexture(const PropertyList& props) {
    m_filepath = props.getString("filepath", "");
    //m_load_texture_successed = Load_BMP_Texture(m_filepath, 3, m_width, m_height);
    m_pixels = stbi_load(m_filepath.c_str(), &m_width, &m_height, &m_channels, 3);
    m_channels = 3;
    if(m_pixels){ m_load_texture_successed = true; }
    else { 
        std::cout << "The image texture did not successfully loaded!" << std::endl;
        m_load_texture_successed = false; 
    }
}


template <>
std::string ImageTexture<Color3f>::toString() const {
    return tfm::format(
        "ImageTexture[\n"
        "  Load file = %s,\n"
        "]",
        m_filepath
    );
}

NORI_REGISTER_TEMPLATED_CLASS(ImageTexture, Color3f, "imagetexture_color3f")
NORI_NAMESPACE_END