#include <unordered_map>
#include <unordered_set>
#include <variant>

class GLFWwindow;

enum class EventState {
    kInit = 0,
    kBuild = 0,
    kSimulate = 0,
};

using KeyHandler = std::function<void(GLFWwindow*, int)>;
using MouseHandler = std::function<void(GLFWwindow*)>;
using MouseMoveHandler = std::function<void(GLFWwindow*, double, double)>;

//
// #############################################################################
//

namespace detail {
enum class Type : uint8_t {
    kSingle = 0,
    kDouble = 1,
    kHold = 2,
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
    Type type = Type::kSingle;
    int modifiers = 0;
};

template <class... Ts>
struct Overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts>
Overloaded(Ts...) -> Overloaded<Ts...>;
}

//
// #############################################################################
//

class EventBuilder {
public:
    // Modifiers
    EventBuilder& control()
    {
        data_.modifiers |= GLFW_MOD_CONTROL;
        return *this;
    }
    EventBuilder& shift()
    {
        data_.modifiers |= GLFW_MOD_SHIFT;
        return *this;
    }

    // Type
    EventBuilder& twice()
    {
        data_.type = detail::Type::kDouble;
        return *this;
    }
    EventBuilder& hold()
    {
        data_.type = detail::Type::kHold;
        return *this;
    }

    // Events
    EventBuilder& key(int key, KeyHandler handler)
    {
        data_.event = detail::KeyEvent { key, std::move(handler) };
        return *this;
    }
    EventBuilder& move(MouseMoveHandler handler)
    {
        data_.event = detail::MouseMoveEvent { std::move(handler) };
        return *this;
    }
    EventBuilder& left_click(MouseHandler handler)
    {
        data_.event = detail::LeftButtonMouseEvent { std::move(handler) };
        return *this;
    }
    EventBuilder& right_click(MouseHandler handler)
    {
        data_.event = detail::RightButtonMouseEvent { std::move(handler) };
        return *this;
    }

    // Get the finished result
    const detail::EventWithMetadata& get_event_with_metadata() const
    {
        return data_;
    }

private:
    detail::EventWithMetadata data_;
};

//
// #############################################################################
//

class EventHandler {
public:
    static EventHandler* get()
    {
        return instance_;
    }

public:
    EventHandler()
    {
        if (instance_) {
            throw std::runtime_error("Only one event handler allowed!");
        }
        instance_ = this;
    }
    ~EventHandler()
    {
        instance_ = nullptr;
    }

public:
    void set_state(EventState state)
    {
        state_ = state;
    }
    EventState get_state() const
    {
        return state_;
    }
    void add_handler(EventState state, const EventBuilder& builder)
    {
        handlers_[state].emplace_back(builder.get_event_with_metadata());
    }

public:
    void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
    {
        if (is_holding_ && action == GLFW_RELEASE) {
            is_holding_ = false;
        }

        const bool double_click = (last_mouse_click_ - Clock::now()) <= kDoubleClick;

        for (const auto& handler : handlers_[state_]) {

            // If this is a double click, don't try to do anything with non-double click handlers
            if (double_click && handler.type != detail::Type::kDouble) {
                continue;
            }

            // Make sure all of the modifiers are met
            if ((handler.modifiers xor mods) != 0) {
                continue;
            }

            // Now we can actually dispatch
            switch (button) {
            case GLFW_MOUSE_BUTTON_LEFT:
                if (auto* event = std::get_if<detail::LeftButtonMouseEvent>(&handler.event))
                    event->callback(window);
                break;
            case GLFW_MOUSE_BUTTON_RIGHT:
                if (auto* event = std::get_if<detail::RightButtonMouseEvent>(&handler.event))
                    event->callback(window);
                break;
            default:
                break;
            }
        }

        // Update this here so that next cycle we'll have started the hold, otherwise a single click will be a hold
        if (action == GLFW_PRESS) {
            is_holding_ = true;
            last_mouse_click_ = Clock::now();
        }
    }

    void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
    {
    }

    void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
    {
    }

private:
    EventState state_;
    std::unordered_map<EventState, std::vector<detail::EventWithMetadata>> handlers_;

    static constexpr auto kDoubleClick = std::chrono::milliseconds(5);
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point last_mouse_click_;
    bool is_holding_ = false;

    static EventHandler* instance_;
};

//
// #############################################################################
//

void hookup_mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
}

void hookup_cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
}

void hookup_key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
}
