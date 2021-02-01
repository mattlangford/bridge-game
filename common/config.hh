#pragma once

namespace common {
/// Canvas size
static constexpr size_t kWidth = 1280;
static constexpr size_t kHeight = 720;

/// How big each block is for rendering
static constexpr size_t kPxSize = 20;  // px
/// How big each block is for simulation
static constexpr size_t kBlockSize = 5;  // meters

/// How many blocks there are on the screen
static constexpr size_t kNumWBlocks = kWidth / kPxSize;
static constexpr size_t kNumHBlocks = kHeight / kPxSize;
}  // namespace common
