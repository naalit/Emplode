#include "Emplode.hpp"

using namespace emplode;

class MyObject : public EmplodeType {
private:
    int counter = 0;
    std::string message = "message 1";

    void PrintMessage() {
        std::cout << message << ": " << counter << std::endl;
    }

    std::string GetMessage() {
        std::stringstream ss;
        ss << message << ": " << counter;
        return ss.str();
    }

    void SetMessage(const std::string &message) {
        this->message = message;
        counter = 0;
    }

public:
    static void InitType(emplode::TypeInfo & info) {
      info.AddMemberFunction("PrintMessage", [](MyObject & target) { target.PrintMessage(); return target.GetMessage(); },
                             "Prints the message represented by the MyObject instance.");
      info.AddMemberFunction("IncCounter", [](MyObject & target) { return ++target.counter; },
                             "Prints the message represented by the MyObject instance.");
    }

    void SetupConfig() override {
        LinkVar(counter, "counter", "The counter");
        LinkFuns<std::string>([&]() { return GetMessage(); }, [&](auto &s) { SetMessage(s); }, "message", "The message");
    }
};

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: Emplode <script>" << std::endl;
        return 1;
    }

    Emplode emplode;
    std::vector<MyObject> objects;
    emplode.AddType<MyObject>("MyObject", "Test object",
        [&](const std::string &name) { objects.push_back(MyObject()); return &objects.back(); },
        [&](const EmplodeType & from, EmplodeType & to) { return true; }
    );
    emplode.Load(argv[1]);

    return 0;
}