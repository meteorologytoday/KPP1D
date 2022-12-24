mutable struct Model

    ev :: Env
    st :: State
    co :: Core
    fo :: Forcing

    function Model(
        ev :: Env,
    )

        st = State(ev)
        co = Core(ev)
        fo = Forcing()



        return new(
            ev,
            st,
            co,
            fo,
        )
    end 


end
