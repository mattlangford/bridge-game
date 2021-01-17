#pragma once
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <variant>

class GLFWwindow;

enum class EventState {
    kInit = 0,
    kBuild = 1,
    kSimulate = 2,
};

using KeyHandler = std::function<void(GLFWwindow*, int)>;
using MouseHandler = std::function<void(GLFWwindow*)>;
using MouseMoveHandler = std::function<void(GLFWwindow*, double, double)>;

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
    EventBuilder& any_modifier()
    {
        data_.modifiers = -1;
        return *this;
    }

    // Type
    EventBuilder& twice()
    {
        data_.type = detail::Type::kDoublePress;
        return *this;
    }
    EventBuilder& hold()
    {
        data_.type = detail::Type::kHold;
        return *this;
    }
    EventBuilder& any_type()
    {
        data_.type = detail::Type::kAny;
        return *this;
    }
    EventBuilder& release()
    {
        data_.type = detail::Type::kRelease;
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

protected:
    // Get the finished result
    const detail::EventWithMetadata& get_event_with_metadata()
    {
        if (std::visit([](const auto& event) { return static_cast<bool>(event.callback); }, data_.event) == false) {
            throw std::runtime_error("Using EventBuilder without an event");
        }
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

private:
    inline static EventHandler* instance_;

    struct BuilderDispatch : EventBuilder {
        template <typename... Args>
        BuilderDispatch(Args&&... args)
            : end(std::forward<Args>(args)...)
        {
        }

        ~BuilderDispatch()
        {
            end(get_event_with_metadata());
        }

        std::function<void(const detail::EventWithMetadata&)> end;
    };

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
    BuilderDispatch add()
    {
        return BuilderDispatch { [this](const detail::EventWithMetadata& data) {
            // TODO add some kind of global handler
            handlers_[EventState::kInit].emplace_back(data);
            handlers_[EventState::kBuild].emplace_back(data);
            handlers_[EventState::kSimulate].emplace_back(data);
        } };
    }
    BuilderDispatch add(EventState state)
    {
        return BuilderDispatch { [state, this](const detail::EventWithMetadata& data) {
            handlers_[state].emplace_back(data);
        } };
    }

public:
    void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
    {
        mods_ = mods;
        if (is_holding_ && action == GLFW_RELEASE) {
            is_holding_ = false;
            return;
        }

        const bool double_click = (Clock::now() - last_mouse_click_) <= kDoubleClick;

        for (const auto& handler : handlers_[state_]) {
            if (handler.type != detail::Type::kAny) {
                using Type = detail::Type;
                auto& type = handler.type;

                if ((type == Type::kHold) xor is_holding_)
                    continue;
                if (type == Type::kRelease && action != GLFW_RELEASE)
                    continue;
                if ((type == Type::kPress || type == Type::kDoublePress) && action != GLFW_PRESS)
                    continue;
            }

            // Make sure all of the modifiers are met
            if ((handler.modifiers >= 0) && (handler.modifiers xor mods_) != 0) {
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
        for (const auto& handler : handlers_[state_]) {
            if ((handler.type != detail::Type::kAny) && (handler.type == detail::Type::kHold) xor is_holding_) {
                continue;
            }

            // Make sure all of the modifiers are met
            if ((handler.modifiers >= 0) && (handler.modifiers xor mods_) != 0) {
                continue;
            }

            if (auto* event = std::get_if<detail::MouseMoveEvent>(&handler.event)) {
                event->callback(window, xpos, ypos);
            }
        }
    }

    void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
    {
        mods_ = mods;

        for (const auto& handler : handlers_[state_]) {
            if (handler.type != detail::Type::kAny) {
                using Type = detail::Type;
                auto& type = handler.type;

                if (type == Type::kHold && action != GLFW_REPEAT)
                    continue;
                if (type == Type::kRelease && action != GLFW_RELEASE)
                    continue;
                if ((type == Type::kPress || type == Type::kDoublePress) && action != GLFW_PRESS)
                    continue;
            }

            // Make sure all of the modifiers are met
            if ((handler.modifiers >= 0) && (handler.modifiers xor mods_) != 0) {

                // Control or shift modifiers are fine, if that's the key being pressed. Only skip if the key isn't one
                // of those two
                bool control_key = key == GLFW_KEY_LEFT_CONTROL || key == GLFW_KEY_RIGHT_CONTROL;
                bool shift_key = key == GLFW_KEY_LEFT_SHIFT || key == GLFW_KEY_RIGHT_SHIFT;
                if (((mods_ xor GLFW_MOD_CONTROL) == 0 && !control_key) || ((mods_ xor GLFW_MOD_SHIFT) == 0 && !shift_key)) {
                    continue;
                }
            }

            // Now we can actually dispatch
            if (auto* event = std::get_if<detail::KeyEvent>(&handler.event)) {
                if (event->key == key) {
                    event->callback(window, key);
                }
            }
        }
    }

private:
    EventState state_ = EventState::kInit;

    std::unordered_map<EventState, std::vector<detail::EventWithMetadata>> handlers_;

    static constexpr auto kDoubleClick = std::chrono::milliseconds(200);
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point last_mouse_click_;
    bool is_holding_ = false;

    int mods_ = 0;
};

//
// #############################################################################
//

void route_mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (auto* handler = EventHandler::get()) {
        handler->mouse_button_callback(window, button, action, mods);
    }
}

void route_cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (auto* handler = EventHandler::get()) {
        handler->cursor_position_callback(window, xpos, ypos);
    }
}

void route_key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (auto* handler = EventHandler::get()) {
        handler->key_callback(window, key, scancode, action, mods);
    }
}
