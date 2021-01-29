#include "renderer/events.hh"

#include <GLFW/glfw3.h>

#include <vector>

namespace {
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

std::string event_state_to_string(const EventState& event_state)
{
    switch (event_state) {
    case EventState::kInit:
        return "Init";
    case EventState::kBuild:
        return "Build";
    case EventState::kSimulate:
        return "Simulate";
    default:
        return "unknown";
    }
}

//
// #############################################################################
//

EventBuilder& EventBuilder::control()
{
    data_.modifiers |= GLFW_MOD_CONTROL;
    return *this;
}
EventBuilder& EventBuilder::shift()
{
    data_.modifiers |= GLFW_MOD_SHIFT;
    return *this;
}
EventBuilder& EventBuilder::any_modifier()
{
    data_.modifiers = -1;
    return *this;
}

//
// #############################################################################
//

EventBuilder& EventBuilder::twice()
{
    data_.type = detail::Type::kDoublePress;
    return *this;
}
EventBuilder& EventBuilder::hold()
{
    data_.type = detail::Type::kHold;
    return *this;
}
EventBuilder& EventBuilder::any_type()
{
    data_.type = detail::Type::kAny;
    return *this;
}
EventBuilder& EventBuilder::release()
{
    data_.type = detail::Type::kRelease;
    return *this;
}

//
// #############################################################################
//

EventBuilder& EventBuilder::key(int key, KeyHandler handler)
{
    data_.event = detail::KeyEvent { key, std::move(handler) };
    return *this;
}
EventBuilder& EventBuilder::move(MouseMoveHandler handler)
{
    data_.event = detail::MouseMoveEvent { std::move(handler) };
    return *this;
}
EventBuilder& EventBuilder::left_click(MouseHandler handler)
{
    data_.event = detail::LeftButtonMouseEvent { std::move(handler) };
    return *this;
}
EventBuilder& EventBuilder::right_click(MouseHandler handler)
{
    data_.event = detail::RightButtonMouseEvent { std::move(handler) };
    return *this;
}

//
// #############################################################################
//

const detail::EventWithMetadata& EventBuilder::get_event_with_metadata()
{
    if (std::visit([](const auto& event) { return static_cast<bool>(event.callback); }, data_.event) == false) {
        throw std::runtime_error("Using EventBuilder without an event");
    }
    return data_;
}

//
// #############################################################################
//

EventHandler::EventHandler()
{
    if (instance_) {
        throw std::runtime_error("Only one event handler allowed!");
    }
    instance_ = this;
    set_state(EventState::kInit);
}

//
// #############################################################################
//

EventHandler::~EventHandler()
{
    instance_ = nullptr;
}

//
// #############################################################################
//

void EventHandler::set_state(EventState state)
{
    if (get_state() == state)
        return;
    state_ = state;
    for (auto& callback : state_callbacks_[state_]) {
        callback(state_);
    }
}

//
// #############################################################################
//

EventState EventHandler::get_state() const
{
    return state_;
}

//
// #############################################################################
//

void EventHandler::add_state_callback(StateChange callback)
{
    // TODO add some kind of global handler
    state_callbacks_[EventState::kInit].emplace_back(callback);
    state_callbacks_[EventState::kBuild].emplace_back(callback);
    state_callbacks_[EventState::kSimulate].emplace_back(callback);
}

//
// #############################################################################
//

void EventHandler::add_state_callback(EventState state, StateChange callback)
{
    state_callbacks_[state].emplace_back(std::move(callback));
}

//
// #############################################################################
//

auto EventHandler::add() -> BuilderDispatch
{
    return BuilderDispatch { [this](const detail::EventWithMetadata& data) {
        // TODO add some kind of global handler
        handlers_[EventState::kInit].emplace_back(data);
        handlers_[EventState::kBuild].emplace_back(data);
        handlers_[EventState::kSimulate].emplace_back(data);
    } };
}

//
// #############################################################################
//

auto EventHandler::add(EventState state) -> BuilderDispatch
{
    return BuilderDispatch { [state, this](const detail::EventWithMetadata& data) {
        handlers_[state].emplace_back(data);
    } };
}

//
// #############################################################################
//

void EventHandler::mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    mods_ = mods;
    if (is_holding_ && action == GLFW_RELEASE) {
        is_holding_ = false;
        return;
    }

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

//
// #############################################################################
//

void EventHandler::cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
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

//
// #############################################################################
//

void EventHandler::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
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
