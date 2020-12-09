
#include <cstddef>
#include <cstdio>
#include "jpeglib.h"
#include "jerror.h"

typedef struct {
  struct jpeg_error_mgr pub;  // "public" fields for IJG library
  // jmp_buf setjmp_buffer;      // For handling catastropic errors
} decoder_error_mgr;

#include "test_bytes.h"

#define RLBOX_MEASURE_TRANSITION_TIMES
#define RLBOX_SINGLE_THREADED_INVOCATIONS
#include "rlbox_types.hpp"

#ifdef MOZ_WASM_SANDBOXING_JPEG
    #if defined(MOZ_WASM_SANDBOXING_MPKFULLSAVE) || defined(MOZ_WASM_SANDBOXING_MPKZEROCOST)
        #  include "rlbox_mpk_sandbox.hpp"
        using rlbox_jpeg_sandbox_type = rlbox::rlbox_mpk_sandbox;
    #else
        #  include "rlbox_lucet_sandbox.hpp"
        using rlbox_jpeg_sandbox_type = rlbox::rlbox_lucet_sandbox;
    #endif
#else
// Extra configuration for no-op sandbox
#  define RLBOX_USE_STATIC_CALLS() rlbox_noop_sandbox_lookup_symbol
#  include "rlbox_noop_sandbox.hpp"
using rlbox_jpeg_sandbox_type = rlbox::rlbox_noop_sandbox;
#endif

#include "rlbox.hpp"

using rlbox_sandbox_jpeg = rlbox::rlbox_sandbox<rlbox_jpeg_sandbox_type>;
template <typename T>
using sandbox_callback_jpeg = rlbox::sandbox_callback<T, rlbox_jpeg_sandbox_type>;
template <typename T>
using tainted_jpeg = rlbox::tainted<T, rlbox_jpeg_sandbox_type>;
template <typename T>
using tainted_opaque_jpeg = rlbox::tainted_opaque<T, rlbox_jpeg_sandbox_type>;
template <typename T>
using tainted_volatile_jpeg = rlbox::tainted_volatile<T, rlbox_jpeg_sandbox_type>;
using rlbox::tainted_boolean_hint;

// Struct info needed for rlbox_load_structs_from_library
extern "C" {
#include "jpeglib.h"
}

#include "JpegStructsForRLBox.h"
rlbox_load_structs_from_library(jpeg);

void MOZ_ASSERT(bool cond, const char* msg = "") {
    if (!cond) {
        std::cout << "Error: " << msg << "\n";
        abort();
    }
}

#define MOZ_RELEASE_ASSERT(...) MOZ_ASSERT(__VA_ARGS__)
#define NS_ASSERTION(...) MOZ_ASSERT(__VA_ARGS__)
#define MOZ_ASSERT_UNREACHABLE(msg) MOZ_ASSERT(false, "Unreachable: " msg)

#define MOZ_LOG(...)
#define LOG_SCOPE(...)


/******************************************************************************/

tainted_opaque_jpeg<unsigned char*> transfer_input_bytes(
  rlbox_sandbox_jpeg* mSandbox,
  unsigned char* buffer, size_t size,
  tainted_opaque_jpeg<unsigned char*>& transfer_buffer,
  size_t& transfer_buffer_size,
  bool& used_copy)
{
  if (transfer_buffer_size >= size) {
    used_copy = true;
    return transfer_buffer;
  } else if (transfer_buffer_size != 0) {
    mSandbox->free_in_sandbox(transfer_buffer);
    transfer_buffer_size = 0;
  }

  const bool free_src_on_copy = false;
  auto transferred = rlbox::copy_memory_or_grant_access(*mSandbox, buffer, size, free_src_on_copy, used_copy);
  MOZ_RELEASE_ASSERT(transferred != nullptr);

  if (used_copy) {
    transfer_buffer = transferred.to_opaque();
    transfer_buffer_size = size;
    return transfer_buffer;
  } else {
    return transferred.to_opaque();
  }
}

tainted_opaque_jpeg<unsigned char*> transfer_input_bytes(
rlbox_sandbox_jpeg* mSandbox,
  unsigned char* buffer, size_t size,
  tainted_opaque_jpeg<unsigned char*>& transfer_buffer,
  size_t& transfer_buffer_size)
{

  bool used_copy = false;
  return transfer_input_bytes(mSandbox, buffer, size, transfer_buffer, transfer_buffer_size, used_copy);
}

/******************************************************************************/

tainted_opaque_jpeg<unsigned char*> m_input_transfer_buffer;
size_t m_input_transfer_buffer_size = 0;
size_t curr_input_ptr_idx = 0;
size_t remaining_bytes = sizeof(inputData);

METHODDEF(tainted_jpeg<boolean>)
fill_input_buffer(rlbox_sandbox_jpeg& aSandbox, tainted_jpeg<j_decompress_ptr> jd) {
    tainted_jpeg<jpeg_source_mgr*> src = jd->src;
    auto to_transfer = remaining_bytes > 400? 400 : remaining_bytes;
    remaining_bytes -= to_transfer;

    auto new_buflen = to_transfer;
    auto transferred = transfer_input_bytes(&aSandbox, const_cast<JOCTET *>(&(inputData[curr_input_ptr_idx])), new_buflen,
        m_input_transfer_buffer, m_input_transfer_buffer_size);
    src->next_input_byte = rlbox::from_opaque(transferred);
    src->bytes_in_buffer = (size_t)new_buflen;

    curr_input_ptr_idx += to_transfer;
    return true;
}

METHODDEF(void)
init_source(rlbox_sandbox_jpeg& aSandbox, tainted_jpeg<j_decompress_ptr> jd) {}

METHODDEF(void)
skip_input_data(rlbox_sandbox_jpeg& aSandbox, tainted_jpeg<j_decompress_ptr> jd, tainted_jpeg<long> num_bytes) {
  tainted_jpeg<jpeg_source_mgr*> src = jd->src;

  if ((num_bytes > rlbox::sandbox_static_cast<long>(src->bytes_in_buffer)).unverified_safe_because(
    "Branches either set tainted data or mBytesToSkip which is checked")) {
      abort();
  } else {
    // Simple case. Just advance buffer pointer

    src->bytes_in_buffer -= rlbox::sandbox_static_cast<size_t>(num_bytes);
    src->next_input_byte += num_bytes;
  }
}

METHODDEF(void)
my_error_exit(rlbox_sandbox_jpeg& aSandbox, tainted_jpeg<j_common_ptr> cinfo) {
  tainted_jpeg<decoder_error_mgr*> err = rlbox::sandbox_reinterpret_cast<decoder_error_mgr*>(cinfo->err);

  // Convert error to a browser error code
  bool memError = err->pub.msg_code.unverified_safe_because("Only checking to set an error code") == JERR_OUT_OF_MEMORY;
  if (memError) {
      printf("OOM error\n");
  }

#ifdef DEBUG
  // char buffer[JMSG_LENGTH_MAX];

  // Create the message
  //(*err->pub.format_message)(cinfo, buffer);

  fprintf(stderr, "JPEG decoding error");
#endif

  abort();
}

METHODDEF(void)
term_source(rlbox_sandbox_jpeg& aSandbox, tainted_jpeg<j_decompress_ptr> jd) {}

/******************************************************************************/

static inline constexpr char RLBOX_JPEG_STATE_ASSERTION[] =
    "Tainted data is being inspected only to check the internal state of "
    "libogg structures. This is not a condition that is critical for safety of "
    "the renderer.";

typedef enum {
  JPEG_HEADER,  // Reading JFIF headers
  JPEG_START_DECOMPRESS,
  JPEG_DECOMPRESS_PROGRESSIVE,  // Output progressive pixels
  JPEG_DECOMPRESS_SEQUENTIAL,   // Output sequential pixels
  JPEG_DONE,
  JPEG_SINK_NON_JPEG_TRAILER,  // Some image files have a
                               // non-JPEG trailer
  JPEG_ERROR
} jstate;

jstate mState = JPEG_HEADER;

enum class WriteState : uint8_t {
  NEED_MORE_DATA,  /// The lambda ran out of data.

  FINISHED,  /// The lambda is done writing to the surface; future writes
             /// will fail.

  FAILURE  /// The lambda encountered an error. The caller may recover
           /// if possible and continue to write. (This never indicates
           /// an error in the SurfacePipe machinery itself; it's only
           /// generated by the lambdas.)
};

#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && \
    defined(__ORDER_BIG_ENDIAN__)
#  if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#    define MOZ_LITTLE_ENDIAN() 1
#    define MOZ_BIG_ENDIAN() 0
#  elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#    define MOZ_LITTLE_ENDIAN() 0
#    define MOZ_BIG_ENDIAN() 1
#  else
#    error "Can't handle mixed-endian architectures"
#  endif
#else
#  error "Don't know how to determine endianness"
#endif

// Represents the bit-shifts required to access color channels when the layout
// is viewed as a uint32_t value.
enum class SurfaceFormatBit : uint32_t {
#if MOZ_LITTLE_ENDIAN()
  R8G8B8A8_R = 0,
  R8G8B8A8_G = 8,
  R8G8B8A8_B = 16,
  R8G8B8A8_A = 24,
#elif MOZ_BIG_ENDIAN()
  R8G8B8A8_A = 0,
  R8G8B8A8_B = 8,
  R8G8B8A8_G = 16,
  R8G8B8A8_R = 24,
#else
#  error "bad endianness"
#endif

  // The following values are endian-independent for A8R8G8B8_UINT32.
  A8R8G8B8_UINT32_B = 0,
  A8R8G8B8_UINT32_G = 8,
  A8R8G8B8_UINT32_R = 16,
  A8R8G8B8_UINT32_A = 24,

  // The following values are OS and endian-independent.
  //
  // TODO(aosmond): When everything blocking bug 1581828 has been resolved, we
  // can make this use R8G8B8A8_X for non-Windows platforms.
  OS_R = A8R8G8B8_UINT32_R,
  OS_G = A8R8G8B8_UINT32_G,
  OS_B = A8R8G8B8_UINT32_B,
  OS_A = A8R8G8B8_UINT32_A,
};

inline uint32_t operator<<(uint8_t a, SurfaceFormatBit b) {
  return a << static_cast<uint32_t>(b);
}

inline uint32_t operator>>(uint32_t a, SurfaceFormatBit b) {
  return a >> static_cast<uint32_t>(b);
}

static void cmyk_convert_bgra(uint32_t* aInput, uint32_t* aOutput,
                              int32_t aWidth) {
  uint8_t* input = reinterpret_cast<uint8_t*>(aInput);

  for (int32_t i = 0; i < aWidth; ++i) {
    const uint32_t iC = input[0];
    const uint32_t iM = input[1];
    const uint32_t iY = input[2];
    const uint32_t iK = input[3];

    const uint8_t r = iC * iK / 255;
    const uint8_t g = iM * iK / 255;
    const uint8_t b = iY * iK / 255;

    *aOutput++ = (0xFF << SurfaceFormatBit::OS_A) |
                 (r << SurfaceFormatBit::OS_R) | (g << SurfaceFormatBit::OS_G) |
                 (b << SurfaceFormatBit::OS_B);
    input += 4;
  }
}

uint32_t* mCMSLine = nullptr;
size_t mCMSLineWidth = 0;
tainted_opaque_jpeg<unsigned char*> m_output_transfer_buffer;
size_t m_output_transfer_buffer_size = 0;

tainted_opaque_jpeg<unsigned char**> m_p_output_transfer_buffer;

uint32_t* aPixelBlock;
int32_t aBlockSize;

WriteState OutputScanlines(
    rlbox_sandbox_jpeg* mSandbox,
    tainted_volatile_jpeg<jpeg_decompress_struct>& mInfo,
    tainted_volatile_jpeg<jpeg_source_mgr>& mSourceMgr,
    tainted_volatile_jpeg<decoder_error_mgr>& mErr
) {
  auto l = [&](uint32_t* aPixelBlock, int32_t aBlockSize) {
        JSAMPROW sampleRow = (JSAMPROW)(mCMSLine ? mCMSLine : aPixelBlock);

        bool used_copy = false;
        auto row_size = (mInfo.output_width * mInfo.output_components).UNSAFE_unverified();
        auto output_buffer = transfer_input_bytes(mSandbox, sampleRow, row_size, m_output_transfer_buffer, m_output_transfer_buffer_size, used_copy);
        auto t_output_buffer = rlbox::from_opaque(output_buffer);
        *rlbox::from_opaque(m_p_output_transfer_buffer) = t_output_buffer;

        if (sandbox_invoke(*mSandbox, jpeg_read_scanlines, &mInfo, m_p_output_transfer_buffer, 1).UNSAFE_unverified() != 1) {
          abort();
        }

        if (used_copy) {
          memcpy(sampleRow, t_output_buffer.UNSAFE_unverified(), row_size);
        }

        switch (mInfo.out_color_space.UNSAFE_unverified()) {
          default:
            // Already outputted directly to aPixelBlock as BGRA.
            MOZ_ASSERT(!mCMSLine);
            break;
          case JCS_GRAYSCALE:
            abort();
            break;
          case JCS_CMYK:
            // Convert from CMYK to BGRA
            MOZ_ASSERT(mCMSLine);
            cmyk_convert_bgra(mCMSLine, aPixelBlock, aBlockSize);
            break;
        }

        return WriteState::FINISHED;
      };

    WriteState result;
    do {
        result = l(aPixelBlock, aBlockSize * 4);
        if (result != WriteState::FINISHED) {
            return result;
        }
    } while (mInfo.output_scanline.UNSAFE_unverified() != mInfo.output_height.UNSAFE_unverified());

    return result;
}

void ReadJPEGData(
    rlbox_sandbox_jpeg* mSandbox,
    tainted_volatile_jpeg<jpeg_decompress_struct>& mInfo,
    tainted_volatile_jpeg<jpeg_source_mgr>& mSourceMgr,
    tainted_volatile_jpeg<decoder_error_mgr>& mErr
) {
  MOZ_LOG(sJPEGLog, LogLevel::Debug,
          ("[this=%p] nsJPEGDecoder::Write -- processing JPEG data\n", this));

  switch (mState) {
    case JPEG_HEADER: {
      LOG_SCOPE((mozilla::LogModule*)sJPEGLog,
                "nsJPEGDecoder::Write -- entering JPEG_HEADER"
                " case");

      auto status = sandbox_invoke(*mSandbox, jpeg_read_header, &mInfo, TRUE);
      // Step 3: read file parameters with jpeg_read_header()
      if (status.unverified_safe_because(RLBOX_JPEG_STATE_ASSERTION) == JPEG_SUSPENDED) {
        abort();
      }

      mSandbox->clear_transition_times();

      // We're doing a full decode.
      switch (mInfo.jpeg_color_space.UNSAFE_unverified()) {
        case JCS_GRAYSCALE:
        case JCS_RGB:
        case JCS_YCbCr:
          // By default, we will output directly to BGRA. If we need to apply
          // special color transforms, this may change.
        //   switch (SurfaceFormat::OS_RGBX) {
        //     case SurfaceFormat::B8G8R8X8:
        //       mInfo.out_color_space = JCS_EXT_BGRX;
        //       break;
        //     case SurfaceFormat::X8R8G8B8:
              mInfo.out_color_space = JCS_EXT_XRGB;
        //       break;
        //     case SurfaceFormat::R8G8B8X8:
        //       mInfo.out_color_space = JCS_EXT_RGBX;
        //       break;
        //     default:
        //       mState = JPEG_ERROR;
        //       return Transition::TerminateFailure();
        //   }
          break;
        case JCS_CMYK:
        case JCS_YCCK:
          // libjpeg can convert from YCCK to CMYK, but not to XRGB.
          mInfo.out_color_space = JCS_CMYK;
          break;
        default:
          mState = JPEG_ERROR;
          MOZ_LOG(sJPEGDecoderAccountingLog, LogLevel::Debug,
                  ("} (unknown colorspace (3))"));
          abort();
      }

      // We don't want to use the pipe buffers directly because we don't want
      // any reads on non-BGRA formatted data.
      if (mInfo.out_color_space.UNSAFE_unverified() == JCS_GRAYSCALE ||
          mInfo.out_color_space.UNSAFE_unverified() == JCS_CMYK) {
        mCMSLineWidth = mInfo.image_width.UNSAFE_unverified();
        mCMSLine = new (std::nothrow) uint32_t[mCMSLineWidth];
        if (!mCMSLine) {
          mState = JPEG_ERROR;
          MOZ_LOG(sJPEGDecoderAccountingLog, LogLevel::Debug,
                  ("} (could allocate buffer for color conversion)"));
          abort();
        }
      }

      // Don't allocate a giant and superfluous memory buffer
      // when not doing a progressive decode.
      mInfo.buffered_image = sandbox_invoke(*mSandbox, jpeg_has_multiple_scans, &mInfo);

      /* Used to set up image size so arrays can be allocated */
      sandbox_invoke(*mSandbox, jpeg_calc_output_dimensions, &mInfo);

      MOZ_LOG(sJPEGDecoderAccountingLog, LogLevel::Debug,
              ("        JPEGDecoderAccounting: nsJPEGDecoder::"
               "Write -- created image frame with %ux%u pixels",
               mInfo.image_width.UNSAFE_unverified(), mInfo.image_height.UNSAFE_unverified()));

      mState = JPEG_START_DECOMPRESS;
      [[fallthrough]];  // to start decompressing.
    }

    case JPEG_START_DECOMPRESS: {
      LOG_SCOPE((mozilla::LogModule*)sJPEGLog,
                "nsJPEGDecoder::Write -- entering"
                " JPEG_START_DECOMPRESS case");
      // Step 4: set parameters for decompression

      // FIXME -- Should reset dct_method and dither mode
      // for final pass of progressive JPEG

      mInfo.dct_method = JDCT_ISLOW;
      mInfo.dither_mode = JDITHER_FS;
      mInfo.do_fancy_upsampling = TRUE;
      mInfo.enable_2pass_quant = FALSE;
      mInfo.do_block_smoothing = TRUE;

      // Step 5: Start decompressor
      if (sandbox_invoke(*mSandbox, jpeg_start_decompress, &mInfo).unverified_safe_because(RLBOX_JPEG_STATE_ASSERTION) == FALSE) {
        MOZ_LOG(sJPEGDecoderAccountingLog, LogLevel::Debug,
                ("} (I/O suspension after jpeg_start_decompress())"));
        abort();
      }

      // If this is a progressive JPEG ...
      mState = mInfo.buffered_image.UNSAFE_unverified() ? JPEG_DECOMPRESS_PROGRESSIVE
                                    : JPEG_DECOMPRESS_SEQUENTIAL;
      [[fallthrough]];  // to decompress sequential JPEG.
    }

    case JPEG_DECOMPRESS_SEQUENTIAL: {
      if (mState == JPEG_DECOMPRESS_SEQUENTIAL) {
        LOG_SCOPE((mozilla::LogModule*)sJPEGLog,
                  "nsJPEGDecoder::Write -- "
                  "JPEG_DECOMPRESS_SEQUENTIAL case");

        switch (OutputScanlines(mSandbox, mInfo, mSourceMgr, mErr)) {
          case WriteState::NEED_MORE_DATA:
            MOZ_LOG(
                sJPEGDecoderAccountingLog, LogLevel::Debug,
                ("} (I/O suspension after OutputScanlines() - SEQUENTIAL)"));
            abort();
          case WriteState::FINISHED:
            NS_ASSERTION(mInfo.output_scanline.UNSAFE_unverified() == mInfo.output_height.UNSAFE_unverified(),
                         "We didn't process all of the data!");
            mState = JPEG_DONE;
            break;
          case WriteState::FAILURE:
            mState = JPEG_ERROR;
            MOZ_LOG(sJPEGDecoderAccountingLog, LogLevel::Debug,
                    ("} (Error in pipeline from OutputScalines())"));
            abort();
        }
      }
      [[fallthrough]];  // to decompress progressive JPEG.
    }

    case JPEG_DECOMPRESS_PROGRESSIVE: {
      if (mState == JPEG_DECOMPRESS_PROGRESSIVE) {
        LOG_SCOPE((mozilla::LogModule*)sJPEGLog,
                  "nsJPEGDecoder::Write -- JPEG_DECOMPRESS_PROGRESSIVE case");

        int status;
        do {
          status = sandbox_invoke(*mSandbox, jpeg_consume_input, &mInfo).UNSAFE_unverified();
        } while ((status != JPEG_SUSPENDED) && (status != JPEG_REACHED_EOI));

        while (mState != JPEG_DONE) {
          if (mInfo.output_scanline.UNSAFE_unverified() == 0) {
            int scan = mInfo.input_scan_number.UNSAFE_unverified();

            // if we haven't displayed anything yet (output_scan_number==0)
            // and we have enough data for a complete scan, force output
            // of the last full scan
            if ((mInfo.output_scan_number.UNSAFE_unverified() == 0) && (scan > 1) &&
                (status != JPEG_REACHED_EOI))
              scan--;

            if (!sandbox_invoke(*mSandbox, jpeg_start_output, &mInfo, scan).unverified_safe_because(RLBOX_JPEG_STATE_ASSERTION)) {
              MOZ_LOG(sJPEGDecoderAccountingLog, LogLevel::Debug,
                      ("} (I/O suspension after jpeg_start_output() -"
                       " PROGRESSIVE)"));
              abort();
            }
          }

          if (mInfo.output_scanline.UNSAFE_unverified() == 0xffffff) {
            mInfo.output_scanline = 0;
          }

          switch (OutputScanlines(mSandbox, mInfo, mSourceMgr, mErr)) {
            case WriteState::NEED_MORE_DATA:
              if (mInfo.output_scanline.UNSAFE_unverified() == 0) {
                // didn't manage to read any lines - flag so we don't call
                // jpeg_start_output() multiple times for the same scan
                mInfo.output_scanline = 0xffffff;
              }
              MOZ_LOG(sJPEGDecoderAccountingLog, LogLevel::Debug,
                      ("} (I/O suspension after OutputScanlines() - "
                       "PROGRESSIVE)"));
              abort();
            case WriteState::FINISHED:
              NS_ASSERTION((mInfo.output_scanline == mInfo.output_height).unverified_safe_because(RLBOX_JPEG_STATE_ASSERTION),
                           "We didn't process all of the data!");

              if (!sandbox_invoke(*mSandbox, jpeg_finish_output, &mInfo).unverified_safe_because(RLBOX_JPEG_STATE_ASSERTION)) {
                MOZ_LOG(sJPEGDecoderAccountingLog, LogLevel::Debug,
                        ("} (I/O suspension after jpeg_finish_output() -"
                         " PROGRESSIVE)"));
                abort();
              }

              if (sandbox_invoke(*mSandbox, jpeg_input_complete, &mInfo).UNSAFE_unverified() &&
                  (mInfo.input_scan_number.UNSAFE_unverified() == mInfo.output_scan_number.UNSAFE_unverified())) {
                mState = JPEG_DONE;
              } else {
                mInfo.output_scanline = 0;
              }
              break;
            case WriteState::FAILURE:
              mState = JPEG_ERROR;
              MOZ_LOG(sJPEGDecoderAccountingLog, LogLevel::Debug,
                      ("} (Error in pipeline from OutputScalines())"));
              abort();
          }
        }
      }
      [[fallthrough]];  // to finish decompressing.
    }

    case JPEG_DONE: {
      LOG_SCOPE((mozilla::LogModule*)sJPEGLog,
                "nsJPEGDecoder::ProcessData -- entering"
                " JPEG_DONE case");

      // Step 7: Finish decompression

      if (sandbox_invoke(*mSandbox, jpeg_finish_decompress, &mInfo).unverified_safe_because(RLBOX_JPEG_STATE_ASSERTION) == FALSE) {
        MOZ_LOG(sJPEGDecoderAccountingLog, LogLevel::Debug,
                ("} (I/O suspension after jpeg_finish_decompress() - DONE)"));
        abort();
      }

      // Make sure we don't feed any more data to libjpeg-turbo.
      mState = JPEG_SINK_NON_JPEG_TRAILER;

      // We're done.
      return;
    }
    case JPEG_SINK_NON_JPEG_TRAILER:
      MOZ_LOG(sJPEGLog, LogLevel::Debug,
              ("[this=%p] nsJPEGDecoder::ProcessData -- entering"
               " JPEG_SINK_NON_JPEG_TRAILER case\n",
               this));

      MOZ_ASSERT_UNREACHABLE(
          "Should stop getting data after entering state "
          "JPEG_SINK_NON_JPEG_TRAILER");

      return;

    case JPEG_ERROR:
      MOZ_ASSERT_UNREACHABLE(
          "Should stop getting data after entering state "
          "JPEG_ERROR");

      abort();
  }

  MOZ_ASSERT_UNREACHABLE("Escaped the JPEG decoder state machine");
  abort();
}

int main(int argc, char const *argv[])
{
    rlbox_sandbox_jpeg sandbox;
    sandbox_callback_jpeg<void(*)(jpeg_decompress_struct *)> init_source_cb;
    sandbox_callback_jpeg<void(*)(j_decompress_ptr)> term_source_cb;
    sandbox_callback_jpeg<void(*)(j_decompress_ptr, long)> skip_input_data_cb;
    sandbox_callback_jpeg<boolean(*)(j_decompress_ptr)> fill_input_buffer_cb;
    sandbox_callback_jpeg<void(*)(j_common_ptr)> my_error_exit_cb;

    #ifdef MOZ_WASM_SANDBOXING_JPEG
        #if defined(MOZ_WASM_SANDBOXING_MPKFULLSAVE) || defined(MOZ_WASM_SANDBOXING_MPKZEROCOST)
        sandbox.create_sandbox(mozilla::ipc::GetSandboxedJpegPath().get());
        #else
        // Firefox preloads the library externally to ensure we won't be stopped
        // by the content sandbox
        const bool external_loads_exist = true;
        // See Bug 1606981: In some environments allowing stdio in the wasm sandbox
        // fails as the I/O redirection involves querying meta-data of file
        // descriptors. This querying fails in some environments.
        const bool allow_stdio = false;
        sandbox.create_sandbox(mozilla::ipc::GetSandboxedJpegPath().get(),
                                external_loads_exist, allow_stdio);
        #endif
    #else
        sandbox.create_sandbox();
    #endif
    init_source_cb = sandbox.register_callback(init_source);
    term_source_cb = sandbox.register_callback(term_source);
    skip_input_data_cb = sandbox.register_callback(skip_input_data);
    fill_input_buffer_cb = sandbox.register_callback(fill_input_buffer);
    my_error_exit_cb = sandbox.register_callback(my_error_exit);

    tainted_jpeg<jpeg_decompress_struct*> p_mInfo = sandbox.malloc_in_sandbox<jpeg_decompress_struct>();
    tainted_jpeg<jpeg_source_mgr*> p_mSourceMgr = sandbox.malloc_in_sandbox<jpeg_source_mgr>();
    tainted_jpeg<decoder_error_mgr*> p_mErr = sandbox.malloc_in_sandbox<decoder_error_mgr>();

    auto& mInfo = *p_mInfo;
    auto& mSourceMgr = *p_mSourceMgr;
    auto& mErr = *p_mErr;
    rlbox::memset(sandbox, &mInfo, 0, sizeof(mInfo));
    rlbox::memset(sandbox, &mSourceMgr, 0, sizeof(mSourceMgr));

    auto p_output_transfer_buffer = sandbox.malloc_in_sandbox<unsigned char*>();
    m_p_output_transfer_buffer = p_output_transfer_buffer.to_opaque();

    aBlockSize = 4 * 1024 * 1024;
    aPixelBlock = new uint32_t[aBlockSize];

    mInfo.err = sandbox_invoke(sandbox, jpeg_std_error, &mErr.pub);
    mErr.pub.error_exit = my_error_exit_cb;

    // Step 1: allocate and initialize JPEG decompression object
    sandbox_invoke(sandbox, jpeg_CreateDecompress, &mInfo, JPEG_LIB_VERSION, sizeof(mInfo));
    // Set the source manager
    mInfo.src = &mSourceMgr;

    // Step 2: specify data source (eg, a file)

    // Setup callback functions.
    mSourceMgr.init_source = init_source_cb;
    mSourceMgr.fill_input_buffer = fill_input_buffer_cb;
    mSourceMgr.skip_input_data = skip_input_data_cb;
    mSourceMgr.resync_to_restart = sandbox.get_sandbox_function_address(jpeg_resync_to_restart);
    mSourceMgr.term_source = term_source_cb;

    // Record app markers for ICC data
    for (uint32_t m = 0; m < 16; m++) {
        sandbox_invoke(sandbox, jpeg_save_markers, &mInfo, JPEG_APP0 + m, 0xFFFF);
    }

    ReadJPEGData(&sandbox, *p_mInfo, *p_mSourceMgr, *p_mErr);

    mInfo.src = nullptr;
    sandbox_invoke(sandbox, jpeg_destroy_decompress, &mInfo);

    // Make sure the compiler doesn't get too clever. Do something with outputs
    uint64_t out = 0;
    for (auto i = 0; i < mCMSLineWidth; i++) {
        out += mCMSLine[i];
    }

    for (auto i = 0; i < aBlockSize; i++) {
        out += aPixelBlock[i];
    }

    std::cout << "Byte sum: " << out << "\n";

    delete[] aPixelBlock;
    delete[] mCMSLine;
    sandbox.free_in_sandbox(p_mInfo);
    sandbox.free_in_sandbox(p_mSourceMgr);
    sandbox.free_in_sandbox(p_mErr);
    sandbox.free_in_sandbox(p_output_transfer_buffer);
    init_source_cb.unregister();
    term_source_cb.unregister();
    skip_input_data_cb.unregister();
    fill_input_buffer_cb.unregister();
    my_error_exit_cb.unregister();
    sandbox.destroy_sandbox();
    return 0;
}
