#pragma once
#include <GLFW/glfw3.h>

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <variant>

struct GLFWwindow;

enum class EventState {
    kInit = 0,
    kBuild = 1,
    kSimulate = 2,
};

std::string event_state_to_string(const EventState &event_state);

using KeyHandler = std::function<void(GLFWwindow *, int)>;
using MouseHandler = std::function<void(GLFWwindow *)>;
using MouseMoveHandler = std::function<void(GLFWwindow *, double, double)>;
using StateChange = std::function<void(EventState)>;

//
// #############################################################################
//

namespace detail {
enum class Type : uint8_t {
    kAny = 0,
    kPress = 1,
    kDoublePress = 2,
    kHold = 3,
    kRelease = 4,
};

struct KeyEvent {
    int key;
    KeyHandler callback;
};
struct MouseMoveEvent {
    MouseMoveHandler callback;
};
struct LeftButtonMouseEvent {
    MouseHandler callback;
};
struct RightButtonMouseEvent {
    MouseHandler callback;
};
using Event = std::variant<KeyEvent, MouseMoveEvent, LeftButtonMouseEvent, RightButtonMouseEvent>;

struct EventWithMetadata {
    Event event;
    Type type = Type::kPress;
    int modifiers = 0;
};
}  // namespace detail

//
// #############################################################################
//

class EventBuilder {
   public:
    // Modifiers
    EventBuilder &control();
    EventBuilder &shift();
    EventBuilder &any_modifier();

    // Type
    EventBuilder &twice();
    EventBuilder &hold();
    EventBuilder &any_type();
    EventBuilder &release();

    // Events
    EventBuilder &key(int key, KeyHandler handler);
    EventBuilder &move(MouseMoveHandler handler);
    EventBuilder &left_click(MouseHandler handler);
    EventBuilder &right_click(MouseHandler handler);

   protected:
    // Get the finished result
    const detail::EventWithMetadata &get_event_with_metadata();

   private:
    detail::EventWithMetadata data_;
};

//
// #############################################################################
//

class EventHandler {
   public:
    inline static EventHandler *get() { return instance_; }

   private:
    inline static EventHandler *instance_;

    struct BuilderDispatch : EventBuilder {
        template <typename... Args>
        BuilderDispatch(Args &&...args) : end(std::forward<Args>(args)...) {}

        ~BuilderDispatch() { end(get_event_with_metadata()); }

        std::function<void(const detail::EventWithMetadata &)> end;
    };

   public:
    EventHandler();
    ~EventHandler();

   public:
    void set_state(EventState state);
    EventState get_state() const;

    void add_state_callback(StateChange callback);
    void add_state_callback(EventState state, StateChange callback);

    BuilderDispatch add();
    BuilderDispatch add(EventState state);

   public:
    void mouse_button_callback(GLFWwindow *window, int button, int action, int mods);

    void cursor_position_callback(GLFWwindow *window, double xpos, double ypos);

    void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods);

   private:
    EventState state_ = EventState::kInit;

    std::unordered_map<EventState, std::vector<detail::EventWithMetadata>> handlers_;
    std::unordered_map<EventState, std::vector<StateChange>> state_callbacks_;

    static constexpr auto kDoubleClick = std::chrono::milliseconds(200);
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point last_mouse_click_;
    bool is_holding_ = false;

    int mods_ = 0;
};

//
// #############################################################################
//

void route_mouse_button_callback(GLFWwindow *window, int button, int action, int mods);

void route_cursor_position_callback(GLFWwindow *window, double xpos, double ypos);

void route_key_callback(GLFWwindow *window, int key, int scancode, int action, int mods);
